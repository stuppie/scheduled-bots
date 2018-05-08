"""
This is a generic OWL -> Wikidata importer
Requires a obographs JSON file made from an owl file

"""

import traceback
from functools import partial
from itertools import repeat
from time import strftime, gmtime

import multiprocessing
from tqdm import tqdm
from wikicurie.wikicurie import CurieUtil, default_curie_map
import json
from datetime import datetime

from scheduled_bots import utils
from wikidataintegrator import wdi_core, wdi_property_store, wdi_helpers, wdi_login
from wikidataintegrator.ref_handlers import update_retrieved_if_new
update_retrieved_if_new_P11 = partial(update_retrieved_if_new, retrieved_pid="P11")

from wikidataintegrator.ref_handlers import update_retrieved_if_new_multiple_refs
from wikidataintegrator.wdi_helpers import WikibaseHelper

uri_to_curie = lambda s: s.split("/")[-1].replace("_", ":")
curie_map = default_curie_map.copy()
cu = CurieUtil(curie_map)


class Node:
    def __init__(self, json_node, graph):
        self.id_uri = json_node['id']  # e.g. http://purl.obolibrary.org/obo/DOID_8718
        self.id_curie = uri_to_curie(self.id_uri)
        self.label = json_node.get("lbl")
        self.type = json_node.get("type")
        self.graph = graph
        self.helper = graph.helper
        self.mediawiki_api_url = graph.mediawiki_api_url
        self.sparql_endpoint_url = graph.sparql_endpoint_url
        self.ref_handler = graph.ref_handler

        self.qid = None
        self.xrefs = set()
        self.item = None

        meta = json_node.get("meta", dict())

        if 'xrefs' in meta:
            self.xrefs = set([x['val'] for x in meta['xrefs']])

        namespaces = [x['val'] for x in meta['basicPropertyValues'] if
                      x['pred'] == 'http://www.geneontology.org/formats/oboInOwl#hasOBONamespace'] if \
            'basicPropertyValues' in meta else []
        if namespaces:
            self.namespace = namespaces[0]
        else:
            self.namespace = graph.default_namespace

        self.descr = meta.get('definition', dict()).get('val')
        self.deprecated = meta.get('deprecated', False)
        self.synonyms = set(x['val'] for x in meta.get('synonyms', list()))
        # filter out the label from the synonyms
        self.synonyms.discard(self.label)

        self.qid = None

    def create_ref_statement(self):
        assert self.graph.release_qid, "create the release first (on the graph class)"

        primary_ext_id_pid, primary_ext_id = cu.parse_curie(self.id_curie)
        primary_ext_id_pid = self.helper.get_pid(primary_ext_id_pid)

        stated_in = wdi_core.WDItemID(value=self.graph.release_qid, prop_nr=self.helper.get_pid('P248'),
                                      is_reference=True)
        ref_extid = wdi_core.WDExternalID(value=primary_ext_id, prop_nr=primary_ext_id_pid, is_reference=True)
        ref_retrieved = wdi_core.WDTime(strftime("+%Y-%m-%dT00:00:00Z", gmtime()),
                                        prop_nr=self.helper.get_pid('P813'), is_reference=True)
        reference = [stated_in, ref_retrieved, ref_extid]
        return reference

    def set_label(self, wd_item):
        # only setting the label if its currently blank or a new item is being created
        if wd_item.get_label() == "":
            wd_item.set_label(self.label)

    def set_descr(self, wd_item):
        # if the current description is blank and the new description
        # is something else (and not over 250 characters), use it
        current_descr = wd_item.get_description()
        if current_descr.lower() in {""} and self.descr and len(self.descr) < 250:
            wd_item.set_description(utils.clean_description(self.descr))
        elif current_descr.lower() == "":
            wd_item.set_description(self.graph.DEFAULT_DESCRIPTION)

    def set_aliases(self, wd_item):
        if self.synonyms is not None:
            wd_item.set_aliases(aliases=self.synonyms, append=True)

    def create_statements(self):
        ref = self.create_ref_statement()

        primary_ext_id_pid, primary_ext_id = cu.parse_curie(self.id_curie)
        primary_ext_id_pid = self.helper.get_pid(primary_ext_id_pid)

        # kind of hacky way to make sure this ID is unique
        wdi_property_store.wd_properties[primary_ext_id_pid] = {'core_id': True}

        s = [wdi_core.WDExternalID(primary_ext_id, primary_ext_id_pid, references=[ref])]

        for xref in self.xrefs:
            if xref.split(":")[0] not in cu.curie_map:
                continue
            pid, ext_id = cu.parse_curie(xref)
            pid = self.helper.get_pid(pid)
            s.append(wdi_core.WDExternalID(ext_id, pid, references=[ref]))

        return s

    def create(self, login, write=True):
        # create or get qid
        # creates the primary external ID, the xrefs, instance of (if set), checks label, description, and aliases
        # not other properties (i.e. subclass), as these may require items existing that may not exist yet
        s = self.create_statements()

        primary_ext_id_pid, primary_ext_id = cu.parse_curie(self.id_curie)
        primary_ext_id_pid = self.helper.get_pid(primary_ext_id_pid)

        try:
            self.item = wdi_core.WDItemEngine(
                item_name=self.label, data=s, domain="this doesn't do anything",
                append_value=self.graph.APPEND_PROPS,
                fast_run=self.graph.FAST_RUN,
                fast_run_base_filter=self.graph.FAST_RUN_FILTER,
                fast_run_use_refs=True,
                global_ref_mode='CUSTOM',
                ref_handler=self.ref_handler,
                mediawiki_api_url=self.mediawiki_api_url,
                sparql_endpoint_url=self.sparql_endpoint_url
            )
        except Exception as e:
            traceback.print_exc()
            msg = wdi_helpers.format_msg(primary_ext_id, primary_ext_id_pid, None, str(e), msg_type=type(e))
            wdi_core.WDItemEngine.log("ERROR", msg)
            return
        self.set_label(self.item)
        self.set_descr(self.item)
        self.set_aliases(self.item)

        wdi_helpers.try_write(self.item, record_id=primary_ext_id, record_prop=primary_ext_id_pid,
                              login=login, write=write)

        self.qid = self.item.wd_item_id


def create_node(node, login, write):
    # pickleable function that calls create on node
    node.create(login, write)


class Graph:
    # the following MUST be overridden in a subclass !!!
    NAME = None
    QID = None
    DEFAULT_DESCRIPTION = None
    APPEND_PROPS = None
    FAST_RUN = None
    FAST_RUN_FILTER = None
    PRED_PID_MAP = None

    # the following must be overriden if a custom node type is being used
    NODE_CLASS = Node

    def __init__(self, mediawiki_api_url='https://www.wikidata.org/w/api.php',
                 sparql_endpoint_url='https://query.wikidata.org/sparql'):
        assert self.NAME, "Must initialize subclass"
        self.version = None
        self.date = None
        self.default_namespace = None
        self.json_graph = None
        self.nodes = None
        self.edges = None

        # str: the QID of the release item. e.g.:
        self.release_qid = None
        # dict[str, dict[str, str]]: use for mapping URIs from external ontologies. e.g. {'UBERON': {'1234': 'Q5453'}}
        self.pid_id_mapper = dict()
        # URIs from this owl file mapping to the node python object
        self.uri_node_map = dict()
        # dict[str, str]: mapping of the hasOBONamespace value to its URI (to be used as 'instance of')
        self.namespace_uri = dict()

        self.mediawiki_api_url = mediawiki_api_url
        self.sparql_endpoint_url = sparql_endpoint_url
        self.helper = WikibaseHelper(sparql_endpoint_url)
        # self.engine = wdi_core.WDItemEngine.wikibase_item_engine_factory(mediawiki_api_url, sparql_endpoint_url)

        # self.ref_handler = partial(update_retrieved_if_new, retrieved_pid=self.helper.get_pid("P813"))
        self.ref_handler = update_retrieved_if_new_P11

    def parse_graph(self, json_path, graph_uri):
        with open(json_path) as f:
            d = json.load(f)
        graphs = {g['id']: g for g in d['graphs']}

        self.json_graph = graphs[graph_uri]
        self.parse_meta()
        self.parse_nodes()
        self.filter_nodes()
        self.parse_edges()
        self.filter_edges()
        self.parse_namespace()
        self.nodes = sorted(self.nodes, key=lambda x: x.id_uri)

    def parse_namespace(self):
        # get all of the namespaces used
        namespaces = set(x.namespace for x in self.nodes)

        # get the label to ID map
        label_uri = {x.label: x.id_uri for x in self.nodes}

        # look up the ID for this namespace
        self.namespace_uri = {namespace: label_uri[namespace] for namespace in namespaces}
        print(self.namespace_uri)

    def parse_meta(self):
        meta = self.json_graph['meta']

        # this is a PURI to the release owl file
        # e.g. : http://purl.obolibrary.org/obo/doid/releases/2018-03-02/doid.owl
        self.version = meta['version']
        # convert the version string to a version string to use for the release item in WD
        self.edition = self.version.rsplit("/", 2)[1]
        # date is for the release item as well
        self.date = datetime.strptime(self.edition, '%Y-%m-%d')

        self.default_namespace = [x['val'] for x in meta['basicPropertyValues'] if
                                  x['pred'] == 'http://www.geneontology.org/formats/oboInOwl#default-namespace'][0]

    def parse_nodes(self):
        self.nodes = [self.NODE_CLASS(json_node, self) for json_node in
                      self.json_graph['nodes']]
        self.uri_node_map = {node.id_uri: node for node in self.nodes}

    def parse_edges(self):
        self.edges = self.json_graph['edges']
        # list of dicts, looks like
        # [{'sub': 'http://purl.obolibrary.org/obo/DOID_820',
        # 'pred': 'http://purl.obolibrary.org/obo/RO_0001025',
        # 'obj': 'http://purl.obolibrary.org/obo/UBERON_0000948'}, ...]

    def filter_edges(self):
        # keep ony the edges where the subject is in the list of filtered nodes
        node_uris = {node.id_uri for node in self.nodes}
        self.edges = [edge for edge in self.edges if edge['sub'] in node_uris]
        # and the pred is in the predicate to PID map "PROP_PID_MAP"
        self.edges = [edge for edge in self.edges if edge['pred'] in self.PRED_PID_MAP]

    def filter_nodes(self):
        self.nodes = [x for x in self.nodes if not x.deprecated and x.type == "CLASS"]

    def create_nodes_par(self, login, write=True):
        pool = multiprocessing.Pool(processes=4)
        create_node_f = partial(create_node, login=login, write=write)
        for _ in tqdm(pool.imap(create_node_f, self.nodes, chunksize=1000), total=len(self.nodes), desc="creating items"):
            pass

    def create_nodes(self, login, write=True):
        create_node_f = partial(create_node, login=login, write=write)
        for _ in tqdm(map(create_node_f, self.nodes), total=len(self.nodes), desc="creating items"):
            pass

    def create_release(self, login):
        r = wdi_helpers.Release('{} release {}'.format(self.NAME, self.edition),
                                'Release of the {}'.format(self.NAME), self.edition,
                                archive_url=self.version, edition_of_wdid=self.QID,
                                pub_date=self.date.date().strftime('+%Y-%m-%dT%H:%M:%SZ'),
                                sparql_endpoint_url=self.sparql_endpoint_url,
                                mediawiki_api_url=self.mediawiki_api_url)
        wd_item_id = r.get_or_create(login)
        if wd_item_id:
            self.release_qid = wd_item_id
        else:
            raise ValueError("unable to create release")

    def create_edges(self, login, write=True):
        # get all the edges for a single subject
        # TODO: this skips edges where the subject is not one of our nodes
        for node in tqdm(self.nodes, desc="creating edges"):
            if not node.qid:
                # todo: log this
                continue
            this_uri = node.id_uri
            this_edges = [edge for edge in self.edges if edge['sub'] == this_uri]
            ss = []
            for edge in this_edges:
                s = self.make_statement_from_edge(edge)
                if s and s.get_value():
                    ss.append(s)

            # set instance of using the namespace
            if node.namespace in self.namespace_uri:
                ref = node.create_ref_statement()
                ss.append(wdi_core.WDItemID(self.uri_node_map[self.namespace_uri[node.namespace]].qid,
                                            self.helper.get_pid('P31'), references=[ref]))
            # print("{}".format([(x.get_value(), x.get_prop_nr()) for x in ss]))
            if ss:
                item = wdi_core.WDItemEngine(
                    wd_item_id=node.qid, data=ss, domain="fake news",
                    append_value=self.APPEND_PROPS,
                    fast_run=self.FAST_RUN,
                    fast_run_base_filter=self.FAST_RUN_FILTER,
                    fast_run_use_refs=True,
                    global_ref_mode='CUSTOM',
                    ref_handler=self.ref_handler,
                    sparql_endpoint_url=self.sparql_endpoint_url,
                    mediawiki_api_url=self.mediawiki_api_url
                )
                this_pid, this_value = cu.parse_curie(uri_to_curie(this_uri))
                this_pid = self.helper.get_pid(this_pid)
                wdi_helpers.try_write(item, record_id=this_value, record_prop=this_pid,
                                      login=login, write=write)

    def make_statement_from_edge(self, edge):
        # we can override this to define a custom statement creator that makes a specific
        # statement depending on the edge or whatever else
        print("edge: {}".format(edge))

        # the predicate has to be defined explicitly
        pred_pid = self.PRED_PID_MAP[edge['pred']]

        # The subject is the item that we have a node for
        subj_node = self.uri_node_map[edge['sub']]

        # the object is a URI either in this node or elsewhere
        obj_qid = self.get_object_qid(edge['obj'])

        if obj_qid:
            return wdi_core.WDItemID(obj_qid, pred_pid, references=[subj_node.create_ref_statement()])


    def get_object_qid(self, edge_obj):
        # object in an edge could be anything. it doesn't have to be a URI that exists within this graph
        # for example, we could be running the DO, and it have an object that is an UBERON class

        # first. check if this URI exists in our graph
        obj_qid = None
        if edge_obj in self.uri_node_map:
            return self.uri_node_map[edge_obj].qid

        # if not, check if the prefix exists in wikidata
        obj_pid, obj_value = cu.parse_curie(uri_to_curie(edge_obj))
        # todo: this will throw an error if not found
        obj_pid = self.helper.get_pid(obj_pid)
        # if this property exists, get all of the values for this property
        if obj_pid not in self.pid_id_mapper:
            self.pid_id_mapper[obj_pid] = wdi_helpers.id_mapper(obj_pid, return_as_set=True,
                                                                prefer_exact_match=True,
                                                                endpoint=self.sparql_endpoint_url)

        # look up by the value
        if obj_value in self.pid_id_mapper[obj_pid]:
            obj_qids = self.pid_id_mapper[obj_pid][obj_value]
            if len(obj_qids) == 1:
                return list(obj_qids)[0]
            else:
                print("oh no: {} {}".format(obj_pid, obj_value))
        else:
            print("oh no: {} {}".format(obj_pid, obj_value))

    def __str__(self):
        return "{} {} #nodes:{} #edges:{}".format(self.default_namespace, self.version, len(self.nodes),
                                                  len(self.edges))
