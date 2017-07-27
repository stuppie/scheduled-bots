import pymysql
import paramiko
import pandas as pd
from sshtunnel import SSHTunnelForwarder
from os.path import expanduser

home = expanduser('~')
mypkey = paramiko.RSAKey.from_private_key_file(home + "/.ssh/id_rsa")

sql_hostname = 'wikidatawiki.labsdb'
sql_main_database = 'wikidatawiki_p'
sql_user = "u16054"
sql_pass = "pnuifBwgDTY3xmVG"
ssh_host = 'tools-login.wmflabs.org'
ssh_user = 'gstupp'
ssh_port = 22
sql_port = 3306

"""
query = '''select rev_id, rev_page, rev_timestamp, page_id, page_namespace, page_title, page_touched FROM revision
               inner join page on revision.rev_page = page.page_id WHERE
               rev_timestamp > DATE_FORMAT(DATE_SUB(NOW(),INTERVAL 1 DAY),'%Y%m%d%H%i%s')
               AND rev_user = 282770;'''
"""

def query_wikidata_mysql(query):
    with SSHTunnelForwarder((ssh_host, ssh_port), ssh_username=ssh_user, ssh_pkey=mypkey,
                            remote_bind_address=(sql_hostname, sql_port)) as tunnel:
        conn = pymysql.connect(host='127.0.0.1', user=sql_user, password=sql_pass, db=sql_main_database,
                               port=tunnel.local_bind_port)

        df = pd.read_sql_query(query, conn)
        conn.close()
    return df