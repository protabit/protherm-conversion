import requests, csv, json, urllib2, random
import xml.etree.ElementTree as ET
from Bio import Entrez

def pdb_query(pdb_id):
    """fetch protein information through a PDB lookup"""
    url = "http://www.rcsb.org/pdb/rest/customReport.csv?pdbids={}&customReportColumns=uniprotAcc,uniprotRecommendedName,taxonomy,sequence&service=wsfile&format=csv".format(pdb_id)
    # should return csv with
    #structureId,chainId,uniprotAcc,uniprotRecommendedName,taxonomy,sequence
    prot_dict={}
    f = urllib2.urlopen(url)
    lines = [l for l in csv.reader(f.readlines(), quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL)]
    if len(lines) > 1:
        try:
            prot_dict['id_pdb_id'] = pdb_id
            prot_dict['id_uniprot_id'] = lines[1][2].strip('"')
            prot_dict['id_name'] = lines[1][3].strip('"')
            prot_dict['id_organism'] = lines[1][4].strip('"')
            seqs=[]
            for line in lines[1:]:
                chain_id = line[1]
                sequence = line[5].strip('\n"X')
                seqs.append(chain_id+":"+sequence)
            if len(seqs) > 1:
                prot_dict['id_sequence'] = "/".join(seqs)
            else:
                prot_dict['id_sequence'] = seqs[0][(seqs[0].find(':')+1):]
        except:
            pass
    return prot_dict

def uniprot_query(uniprotid):
    prot_dict = {}
    req = requests.get('http://www.uniprot.org/uniprot/'+uniprotid+'.xml')
    if req.status_code != 200 or len(uniprotid) < 6: # bad status_code likely wrong id
        return prot_dict
    root=ET.fromstring(req.content)
    prot_dict['id_uniprot_id'] = uniprotid
    try:
        prot_dict['id_name'] = root.find('.//{http://uniprot.org/uniprot}fullName').text
    except:
        pass
    try:
        prot_dict['id_sequence'] = root.find('.//{http://uniprot.org/uniprot}sequence').text.replace('\n','')
    except:
        elements = root.findall('.//{http://uniprot.org/uniprot}sequence')
        for element in elements:
            if element.text is not None:
                prot_dict['id_sequence'] = element.text.replace('\n','')
                break
    try:
        prot_dict['id_organism'] = root.find('.//*{http://uniprot.org/uniprot}organism/{http://uniprot.org/uniprot}name').text
    except:
        pass
    try:
        prot_dict['id_pdb_id'] = ",".join(set([el.attrib['id'] for el in root.findall(".//*{http://uniprot.org/uniprot}dbReference/[@type='PDB']")]))
    except:
        pass
    return prot_dict

def pubmed_bulk_query(idlist):
    chunk = 200
    Entrez.email='connie.wang@protabit.com'
    pubmed_dict={}
    for i in range((len(idlist)+199)/chunk):
        print "fetching chunk ",i,"of ", len(idlist)/chunk
        start=i*chunk
        end=(i+1)*chunk
        ids = idlist[start:end]
        handle=Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode='xml')
        record=Entrez.read(handle)
        article_list = record['PubmedArticle']
        for article in article_list:
            try:
                id=str(article['MedlineCitation']['PMID'])
                title=article['MedlineCitation']['Article']['ArticleTitle']
                journal=article['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
                authorlist=article['MedlineCitation']['Article']['AuthorList']
                authors=";".join([author['LastName']+ " "+author['Initials'] for author in authorlist ])
                pubmed_dict[id]={'pubmedid':id, 'abstract':'', 'title':title, 'journal':journal, 'authors':authors, 'volume':'', 'pubdate':''}
                try:
                    pubmed_dict[id]['abstract']=article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
                except:
                    pass
                try:
                    pubmed_dict[id]['issue']=article['MedlineCitation']['Article']['Journal']['JournalIssue']['Issue']
                except:
                    pass
                try:
                    pubmed_dict[id]['pubdate']=article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
                except:
                    pass
                try:
                    pubmed_dict[id]['volume']=article['MedlineCitation']['Article']['Journal']['JournalIssue']['Volume']
                except:
                    pass
                try:
                    pubmed_dict[id]['pages']=article['MedlineCitation']['Article']['Pagination']['MedlinePgn']
                except:
                    pass
            except:
                print "key details missing"
                pass
    return pubmed_dict

def search_for_pubmedid(pub_dict):
    """ match the reference string into journal, issue, year, pages to find the reference """
    Entrez.email='connie.wang@protabit.com'
    term=pub_dict['journal']+"[ta] AND "+pub_dict['volume']+"[vi] AND "+pub_dict['pages']+"[pg]"
    try:
        handle=Entrez.esearch(db="pubmed", term=term)
        record=Entrez.read(handle)
        if record['Count'] == '1':
            return record['IdList'][0]
        else:
            return None
    except:
        print 'Exception', term
        return None

