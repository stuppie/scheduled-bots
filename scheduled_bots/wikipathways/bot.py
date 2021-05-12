import json
import argparse
import copy
import traceback
from datetime import datetime
import os
from rdflib import Graph
import zipfile
import io
from contextlib import closing
from bs4 import BeautifulSoup, SoupStrainer
from SPARQLWrapper import SPARQLWrapper, JSON
import re
import tqdm
import requests
import sys
# from scheduled
from scheduled_bots import PROPS, ITEMS, get_default_core_props
from wikidataintegrator import wdi_core, wdi_login, wdi_helpers
from wikidataintegrator.ref_handlers import update_retrieved_if_new_multiple_refs
from wikidataintegrator.wdi_helpers import try_write

CACHE_SIZE = 10000
CACHE_TIMEOUT_SEC = 300  # 5 min

print("Logging in...")
if "WDUSER" in os.environ and "WDPASS" in os.environ:
    WDUSER = os.environ['WDUSER']
    WDPASS = os.environ['WDPASS']
else:
    raise ValueError("WDUSER and WDPASS must be specified in local.py or as environment variables")

PROPS = {
    'Wikipathways ID': 'P2410',
    'instance of': 'P31',
    'stated in': 'P248',
    'reference URL': 'P854',
    'Entrez Gene ID ': 'P351',
    'found in taxon': 'P703',
    'PubMed ID': 'P698',
    'curator': 'P1640',
    'retrieved': 'P813',
    'medical condition': 'P1050',
    'author': 'P50',
    'author name string': 'P2093',
    'main subject': 'P921',
    'anatomical location': 'P927'
}

core_props = get_default_core_props()
core_props.update({PROPS['Wikipathways ID']})

ITEMS = {
    'Wikipathways': 'Q7999828',
    'Homo sapiens': 'Q15978631'
}

#core_props = get_default_core_props()
#core_props.update({PROPS['WikiPathways ID']})

__metadata__ = {
    'name': 'PathwayBot',
    'maintainer': 'Andra',
    'tags': ['pathways'],
    'properties': list(PROPS.values())
}

fast_run_base_filter = {'P2410': ''}
fast_run = True

# populate ontology dictionaries
# Disease Ontology
doid_qid_query = wdi_core.WDItemEngine.execute_sparql_query("SELECT * WHERE {?qid wdt:P699 ?doid . }")
doid_qid = dict()
for result in doid_qid_query["results"]["bindings"]:
    doid_qid[result["doid"]["value"]] = result["qid"]["value"].replace("http://www.wikidata.org/entity/", "")

# Pathway Ontology
poid_qid_query = wdi_core.WDItemEngine.execute_sparql_query("SELECT * WHERE {?qid wdt:P7333 ?poid . }")
poid_qid = dict()
for result in poid_qid_query["results"]["bindings"]:
    poid_qid[result["poid"]["value"]] = result["qid"]["value"].replace("http://www.wikidata.org/entity/", "")

def create_reference(pathway_id, retrieved):
    refStatedIn = wdi_core.WDItemID(value=ITEMS['Wikipathways'], prop_nr=PROPS['stated in'], is_reference=True)
    timeStringNow = retrieved.strftime("+%Y-%m-%dT00:00:00Z")
    refRetrieved = wdi_core.WDTime(timeStringNow, prop_nr=PROPS['retrieved'], is_reference=True)
    refWikiPathwaysID = wdi_core.WDString(value=pathway_id, prop_nr=PROPS['Wikipathways ID'], is_reference=True)

    pathway_reference = [refStatedIn, refRetrieved, refWikiPathwaysID]
    return pathway_reference

def panic(pathway_id, msg='', msg_type=''):
    s = wdi_helpers.format_msg(pathway_id, PROPS['Wikipathways ID'], None, msg, msg_type)
    wdi_core.WDItemEngine.log("ERROR", s)
    print(s)
    return None

def main(retrieved, fast_run, write):
    login = wdi_login.WDLogin(WDUSER, WDPASS)
    temp = Graph()
    # TODO add try catch to catch api errors while downloading. Currently when there is API glitch, an UTF error is report, which is only because it gets empty content.
    # TODO add progress indicator of the download

    url = 'http://data.wikipathways.org/current/rdf'
    page = requests.get(url).text
    files = []
    for link in BeautifulSoup(page, "lxml", parse_only=SoupStrainer('a')):
        address = str(link).split("\"")
        if len(address) > 1:
            filename = address[1].replace("./", "/")
            if len(filename) > 1:
                if filename not in files:
                    if filename != "./":
                        files.append(url + filename)
    for file in set(files):
        if "rdf-authors" in file:  # get the most accurate file
            print(file)
            u = requests.get(file)
            with closing(u), zipfile.ZipFile(io.BytesIO(u.content)) as archive:
                for member in archive.infolist():
                    if "_" in str(member.filename):
                        continue
                    print("parsing: " + member.filename)
                    nt_content = archive.read(member)
                    temp.parse(data=nt_content, format="turtle")
            print("size: "+str(len(temp)))

        if "rdf-wp" in file:  # get the most accurate file
            print(file)
            u = requests.get(file)
            with closing(u), zipfile.ZipFile(io.BytesIO(u.content)) as archive:
                for member in archive.infolist():
                    nt_content = archive.read(member)
                    temp.parse(data=nt_content.decode(), format="turtle")
            print("size: "+str(len(temp)))


    wp_query = """prefix dcterm: <http://purl.org/dc/terms/>
            prefix wp: <http://vocabularies.wikipathways.org/wp#>
            SELECT DISTINCT ?wpid WHERE {
              ?s rdf:type <http://vocabularies.wikipathways.org/wp#Pathway> ;
                 dcterm:identifier ?wpid ;
                 ?p <http://vocabularies.wikipathways.org/wp#Curation:AnalysisCollection> ;
                 wp:organism <http://purl.obolibrary.org/obo/NCBITaxon_9606> .
              }"""

    wpids = []
    qres = temp.query(wp_query)
    for row in qres:
        wpids.append(str(row[0]))

    for pathway_id in wpids:
        try:
            run_one(pathway_id, retrieved, fast_run, write, login, temp)
        except Exception as e:
            traceback.print_exc()
            wdi_core.WDItemEngine.log("ERROR", wdi_helpers.format_msg(
                pathway_id, PROPS['Wikipathways ID'], None, str(e), type(e)))

def run_one(pathway_id, retrieved, fast_run, write, login, temp):
    print(pathway_id)
    pathway_reference = create_reference(pathway_id, retrieved)
    prep = dict()

    prep = get_PathwayElements(pathway=pathway_id,datatype="Metabolite", temp=temp, prep=prep)
    prep = get_PathwayElements(pathway=pathway_id,datatype="Protein", temp=temp, prep=prep)
    prep = get_PathwayElements(pathway=pathway_id, datatype="GeneProduct",temp=temp, prep=prep)
    prep = get_PathwayElements(pathway=pathway_id, datatype="Complex",temp=temp, prep=prep)
    # P703 = found in taxon, Q15978631 = "Homo sapiens"
    prep["P703"] = [
        wdi_core.WDItemID(value="Q15978631", prop_nr='P703', references=[copy.deepcopy(pathway_reference)])]

    query = """
            PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
            PREFIX gpml:    <http://vocabularies.wikipathways.org/gpml#>
            PREFIX dcterms: <http://purl.org/dc/terms/>
        SELECT DISTINCT ?pathway ?pwId ?pwLabel ?description
        WHERE {
           VALUES ?pwId {"""
    query += "\"" + pathway_id + "\"^^xsd:string}"
    query += """
           ?pathway a wp:Pathway ;
                    dc:title ?pwLabel ;
                    dcterms:description ?description ;
                    dcterms:identifier ?pwId ;
                    <http://vocabularies.wikipathways.org/wp#isAbout> ?details ;
                    wp:organismName "Homo sapiens"^^xsd:string .
        }"""
    qres3 = temp.query(query)

    for row in qres3:
        #pathway_iri = str(row[0])
        pw_id = str(row[1])
        pw_label = str(row[2])
        description = str(row[3])
        
        ## clean up descriptions
        description = re.sub(r'https?:\/\/.*[\s\r\n]','',description)
        description = description.replace('\n',' ').replace('\r',' ').replace('\'\'\'','').replace('\'\'','').replace('[','').replace(']','')
        description = description.replace('Proteins on this pathway have targeted assays available via the Portal','')
        description = (description[:246]+'...') if len(description) > 246 else description
        description = 'biological pathway in human' if len(description) < 20 else description

        # P31 = instance of
        prep["P31"] = [
            wdi_core.WDItemID(value="Q4915012", prop_nr="P31", references=[copy.deepcopy(pathway_reference)])]

        prep["P1476"] = [wdi_core.WDMonolingualText(value=pw_label, prop_nr="P1476", references=[copy.deepcopy(pathway_reference)])]

        # P2410 = WikiPathways ID
        prep["P2410"] = [wdi_core.WDString(pathway_id, prop_nr='P2410', references=[copy.deepcopy(pathway_reference)])]

        # P2888 = exact match
        prep["P2888"] = [wdi_core.WDUrl("http://identifiers.org/wikipathways/" + pw_id, prop_nr='P2888',
                                        references=[copy.deepcopy(pathway_reference)])]

        query = """
                PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
                PREFIX dcterms: <http://purl.org/dc/terms/>
                select ?pubmed

                WHERE {
                 ?pubmed  a       wp:PublicationReference ;
                        dcterms:isPartOf <"""

        query += str(row[0])
        query += """> .}

                """
        qres4 = temp.query(query)
        p = re.compile('^[0-9]+$')
        for pubmed_result in qres4:
            pmid = str(pubmed_result[0]).replace("http://identifiers.org/pubmed/", "").replace("https://identifiers.org/pubmed/", "")
            print(pmid)
            m = p.match(pmid)

            if not m:
                pmid_qid, _, _ = wdi_helpers.PublicationHelper(pmid, id_type="doi",
                                                               source="crossref").get_or_create(login if write else None)
            else:
                pmid_qid, _, _ = wdi_helpers.PublicationHelper(pmid.replace("PMID:", ""), id_type="pmid",
                                                     source="europepmc").get_or_create(login if write else None)
            if pmid_qid  is None:
                return panic(pathway_id, "not found: {}".format(pmid), "pmid")
            else:
                if 'P2860' not in prep.keys():
                    prep["P2860"] = []
                print(pmid_qid)
                prep['P2860'].append(wdi_core.WDItemID(value=str(pmid_qid), prop_nr='P2860', references=[copy.deepcopy(pathway_reference)]))

        author_query = """
                PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
                SELECT ?author ?authorName ?authorHomepage ?authorQIRI
                WHERE {
                  { <http://identifiers.org/wikipathways/""" + pathway_id + """> dc:creator ?author . }
                  UNION
                  { <https://identifiers.org/wikipathways/""" + pathway_id + """> dc:creator ?author . }
                  ?author a                    foaf:Person ;
                  foaf:name            ?authorName ;
                  foaf:homepage            ?authorHomepage .
                  OPTIONAL { ?author    owl:sameAs     ?authorQIRI . }
                }
                """
        author_query_res = temp.query(author_query)
        prep["P2093"] = []
        prep["P50"] = []

        for row in author_query_res:
            author_name = str(row[1])
            print("author_name")
            print(author_name)
            author_homepage = str(row[2])
            print("author_homepage")
            print(author_homepage)

            # P2093 = author name string
            author_url_qualifier = wdi_core.WDString(value=author_homepage, prop_nr="P2699", is_qualifier=True)
            prep["P2093"].append(wdi_core.WDString(author_name, prop_nr='P2093', qualifiers=[copy.deepcopy(author_url_qualifier)], references=[copy.deepcopy(pathway_reference)]))


            if row[3] != None: # only if row[3] exists (authorQIRI)
                author_iri = str(row[0])
                author_name = str(row[1])
                print("author_name")
                print(author_name)
                author_qiri = str(row[3])
                author_qid = author_qiri.replace("http://www.wikidata.org/entity/", "")
                print("author_qid")
                print(author_qid)
                # P50 = author
                prep["P50"].append(wdi_core.WDItemID(author_qid, prop_nr='P50', references=[copy.deepcopy(pathway_reference)]))

        disease_ontology_query = """
                PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
                PREFIX dcterms: <http://purl.org/dc/terms/>
                SELECT ?diseaseOntologyTerm
                WHERE {
                  ?pathwayRDF wp:diseaseOntologyTag ?diseaseOntologyTerm ;
                    foaf:page ?pathway ;
                    dcterms:identifier \"""" + pathway_id + """\"^^xsd:string . 
                }
                """
        disease_ontology_query_res = temp.query(disease_ontology_query)
        prep["P1050"] = []
        for row in disease_ontology_query_res:
            disease_ontology_iri = str(row[0])
            doid = disease_ontology_iri.replace("http://purl.obolibrary.org/obo/DOID_", "DOID:")
            print("doid")
            print(doid)

            # P1050 = medical condition
            if doid_qid.get(doid) != None:  #skip if qid is missing
                prep["P1050"].append(wdi_core.WDItemID(doid_qid[doid], prop_nr='P1050', references=[copy.deepcopy(pathway_reference)]))
            
        pw_ontology_query = """
                PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
                PREFIX dcterms: <http://purl.org/dc/terms/>
                SELECT ?pwOntologyTerm
                WHERE {
                  ?pathwayRDF wp:pathwayOntologyTag ?pwOntologyTerm ;
                    foaf:page ?pathway ;
                    dcterms:identifier \"""" + pathway_id + """\"^^xsd:string . 
                }
                """
        pw_ontology_query_res = temp.query(pw_ontology_query)
        prep["P921"] = []
        for row in pw_ontology_query_res:
            pw_ontology_iri = str(row[0])
            poid = pw_ontology_iri.replace("http://purl.obolibrary.org/obo/PW_", "PW:")
            print("poid")
            print(poid)

            # P921 = main subject
            if poid_qid.get(poid) != None:  #skip if qid is missing
                prep["P921"].append(wdi_core.WDItemID(poid_qid[poid], prop_nr='P921', references=[copy.deepcopy(pathway_reference)]))

        #TODO: Propose Cell Type Ontology ID as new property, add release item, associate terms with WD items.
        #cell_type_ontology_query = """
        #        PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
        #        PREFIX dcterms: <http://purl.org/dc/terms/>
        #        SELECT ?cellTypeOntologyTerm
        #        WHERE {
        #          ?pathwayRDF wp:cellTypeOntologyTag ?cellTypeOntologyTerm ;
        #            foaf:page ?pathway ;
        #            dcterms:identifier \"""" + pathway_id + """\"^^xsd:string . 
        #        }
        #        """
        #cell_type_ontology_query_res = temp.query(cell_type_ontology_query)
        #prep["P927"] = []
        #for row in cell_type_ontology_query_res:
        #    cell_type_ontology_iri = str(row[0])
        #    ctoid = cell_type_ontology_iri.replace("http://purl.obolibrary.org/obo/CL_", "CL:")
        #    print("ctoid")
        #    print(ctoid)

            # P927 = anatomical location
        #    prep["P927"].append(wdi_core.WDItem(qid[ctoid], prop_nr='P927', references=[copy.deepcopy(pathway_reference)]))

        data2add = []
        for key in prep.keys():
            for statement in prep[key]:
                data2add.append(statement)
                print(statement.prop_nr, statement.value)

        wdPage = wdi_core.WDItemEngine(data=data2add,
                                     fast_run=fast_run,
                                     fast_run_base_filter=fast_run_base_filter,
                                     fast_run_use_refs=True,
                                     ref_handler=update_retrieved_if_new_multiple_refs,
                                     core_props=core_props)

        wdPage.set_label(pw_label, lang="en")
        wdPage.set_description(description, lang="en")

        try_write(wdPage, record_id=pathway_id, record_prop=PROPS['Wikipathways ID'],
                edit_summary="Updated a Wikipathways pathway", login=login, write=write)

def get_PathwayElements(pathway, datatype, temp, prep):
    query = """PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>
           PREFIX rdfs:    <http://www.w3.org/2000/01/rdf-schema#>
           PREFIX dcterms: <http://purl.org/dc/terms/>

           select distinct ?pathway (str(?label) as ?geneProduct) ?id where {
            ?metabolite a wp:"""
    query += datatype
    query += """  ;
                    rdfs:label ?label ;"""
    if datatype == "Metabolite":
        query += "   wp:bdbPubChem ?id ;"
    if datatype == "Protein":
        query += "   wp:bdbWikidata ?id ;"
    if datatype == "GeneProduct":
        query += "   wp:bdbEntrezGene ?id ;"
    if datatype == "Complex":
        query += "   wp:bdbWikidata ?id ;"
    query += """
                    dcterms:isPartOf ?pathway .
            ?pathway a wp:Pathway ;
                   dcterms:identifier
            """
    query += "\"" + pathway + "\"^^xsd:string .}"
    qres2 = temp.query(query)

    ids = []
    if datatype == "Protein":
        for row in qres2:
            ids.append("wd:" + str(row[2]).replace(
                "http://www.wikidata.org/entity/", "")
            )
    elif datatype == "Complex":
        for row in qres2:
            ids.append("wd:" + str(row[2]).replace(
                "http://www.wikidata.org/entity/", "")
            )
    else:
        for row in qres2:
            ids.append("\"" + str(row[2]).replace("http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID", "").replace(
                "http://identifiers.org/ncbigene/", "") + "\"")


    # Check for existence of the ids in wikidata

    wd_query = "SELECT DISTINCT * WHERE {VALUES ?id {"
    wd_query += " ".join(list(set(ids)))
    if datatype == "Metabolite":
        wd_query += "} ?item wdt:P662 ?id . }"
    if datatype == "GeneProduct":
        wd_query += "} ?item wdt:P351 ?id . }"
    if datatype == "Protein":
        wd_query = "SELECT DISTINCT * WHERE { VALUES ?item { "
        wd_query += " ".join(list(set(ids)))
        wd_query += " } ?item wdt:P31 | wdt:P279 wd:Q8054 }"
    if datatype == "Complex":
        wd_query = "SELECT DISTINCT * WHERE { VALUES ?item { "
        wd_query += " ".join(list(set(ids)))
        wd_query += " } ?item wdt:P7718 | wdt:P3937 ?id }"

    if datatype == "Complex":
        print(wd_query)

    results = wdi_core.WDItemEngine.execute_sparql_query(wd_query,)
    if datatype == "Complex":
        print(results)
    for result in results["results"]["bindings"]:
        if "P527" not in prep.keys():
            prep["P527"] = []
        pathway_reference = create_reference(pathway, retrieved)
        prep["P527"].append(
            wdi_core.WDItemID(result["item"]["value"].replace("http://www.wikidata.org/entity/", ""),
                              prop_nr='P527', references=[copy.deepcopy(pathway_reference)]))
    return prep

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run wikipathways bot')
    parser.add_argument('--dummy', help='do not actually do write', action='store_true')
    parser.add_argument('--no-fastrun', action='store_true')
    args = parser.parse_args()
    log_dir = "./logs"
    run_id = datetime.now().strftime('%Y%m%d_%H:%M')
    __metadata__['run_id'] = run_id
    fast_run = False if args.no_fastrun else True
    retrieved = datetime.now()

    log_name = '{}-{}.log'.format(__metadata__['name'], run_id)
    if wdi_core.WDItemEngine.logger is not None:
        wdi_core.WDItemEngine.logger.handles = []
    wdi_core.WDItemEngine.setup_logging(log_dir=log_dir, log_name=log_name, header=json.dumps(__metadata__),
                                        logger_name='wikipathways')

    main(retrieved, fast_run=fast_run, write=not args.dummy)
