import re
import os
import requests     # Pre-installation Required

#Open Files
def openfile(file_url):
    name_urls = []
    with open(file_url, 'r', encoding='utf-8') as f:
        contents = f.readlines()
        for content in contents:
            content = content.replace('\t\n', '').replace('\n', '').split('\t')
            name_urls.append(content)
        return name_urls
#Save Files
def savefile(name, res):
    dir = 'D:/Biocontainers/HTML/%s.html'%name
    if not os.path.exists('D:/Biocontainers/HTML/'):    # File Saving Address
        os.mkdir('D:/Biocontainers/HTML/')
    with open(dir, 'w', encoding='utf-8') as f:
        f.write(res)

#Get README.md & save as html files
def getdate(url, headers, i):
    print(url + '   --------------     ' + str(i))
    try:
        res = requests.get(url=url, headers=headers, timeout=15, allow_redirects=False)
        response = res.text
        if ('doi' in response) or ('PMID' in response) or ('Cite' in response) or ('Citation' in response):
            print('       Keyword Exists' + '-----------------'+ str(i))
            savefile(name, response)
        else:
            print('No keyword Exists.' + '-----------------'+ str(i))
    except requests.exceptions.ConnectionError:
        print("Network Error-----------%d" % i)
        url_again.append(url)
    except requests.Timeout:
        print("Request Timeout-----------%d" % i)
        url_again.append(url)
    except requests.exceptions.TooManyRedirects:
        print("Excessive Redirection.")
    except requests.exceptions.InvalidSchema:
        print("Error in the URL.")

    return url_again

if __name__ == "__main__":
    headers = {
        'Accept': 'text/html',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'zh-CN,zh;q=0.9',
        'Connection': 'keep-alive',
        'Cookie': '_octo=GH1.1.2133912576.1600409893; logged_in=no; _ga=GA1.2.1831073498.1600409904; _gh_sess=hS4oMIaMHraJqN9wO6%2BFR0jF0JyPa%2FhotPo0Ujeq5%2BSDuYXAX0Libfa%2BUO%2Bj41BH9NWg1n057KA6LJao9Rtnv6W5C9Nqn2fgJoPDMPXooa8ZSqCtA1sb2LDTCqaezWGf8pz0%2Bmf01pWVrGrgyobbSVHhDT315UDAGl6NJT%2BueaLb%2Fpp%2FCfJ%2B5tSfZsJG8iE%2F6QpAMOGp8088eAd0KNas7fmFwyHcaJd9zJ%2BKjOagBLqqcf%2FfH0MZpuFQeL2yp2oqJlgSfhdmbdrJ83C%2ByvO0Kw%3D%3D--A3%2BF80BsJEuy4LO3--U5vLDuNGDSwCKWHRI9b8GA%3D%3D; tz=Asia%2FShanghai',
        'Host': 'github.com',
        'Referer': 'https://github.com/samtools/samtools',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.105 Safari/537.36',
        'X-Requested-With': 'XMLHttpRequest',
    }
    file_url = 'D:/Biocontainers/missing_dois.tsv'       # Local Address of Yaml Files
    i = 0  # Used for Counting
    url_again = []  # Store Failed Links
    name_urls = openfile(file_url)
    print("%d in total"%len(name_urls))
    for name_url in name_urls:
        i += 1
        name = name_url[0]
        url = name_url[1]
        getdate(url,headers,i)
    if url_again != None:
        i = 0
        print("Try again")
        url_again_len = len(url_again)
        print(url_again_len)
        for url in url_again:
            i += 1
            if i <= url_again_len:
                getdate(url, headers, i)
            else:
                break
