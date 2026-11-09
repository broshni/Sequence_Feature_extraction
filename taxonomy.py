#NCBI Taxonomy lookups shared by the steps that need to know which organism a name refers to.
#Nothing here is specific to one virus: names and IDs are resolved against NCBI, so any virus a
#user searches for works the same way.
import re
from functools import lru_cache

from Bio import Entrez


def configure(email, api_key=None):
    """NCBI requires an email address on every Entrez request."""
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key


def _species_of(entry):
    """The species-level taxid for a Taxonomy entry, or None when it has no species rank."""
    if entry.get('Rank') == 'species':
        return str(entry['TaxId'])
    for ancestor in entry.get('LineageEx', []):
        if ancestor.get('Rank') == 'species':
            return str(ancestor['TaxId'])
    return None


def _read_entries(taxids):
    handle = Entrez.efetch(db='taxonomy', id=','.join(str(t) for t in taxids), retmode='xml')
    entries = Entrez.read(handle)
    handle.close()
    return entries


def describe(entry):
    lineage = [str(entry['TaxId'])] + [str(a['TaxId']) for a in entry.get('LineageEx', [])]
    return {
        'taxid': str(entry['TaxId']),
        'name': str(entry['ScientificName']),
        'rank': str(entry.get('Rank', '')),
        'species_taxid': _species_of(entry),
        'lineage': lineage,
    }


def search(name, retmax=5):
    """Taxonomy IDs NCBI matches to a name, best match first. A strain-style name NCBI has no
    entry for falls back to the name before its first parenthesis."""
    def run(term):
        handle = Entrez.esearch(db='taxonomy', term=f'"{term}"[All Names]', retmax=retmax)
        found = [str(i) for i in Entrez.read(handle)['IdList']]
        handle.close()
        return found
    found = run(name)
    if not found and '(' in name:
        trimmed = re.sub(r'\s*\(.*$', '', name).strip()
        found = run(trimmed) if trimmed else []
    return found


@lru_cache(maxsize=None)
def lookup(query):
    """Resolve a taxonomy ID or an organism name to
    {taxid, name, rank, species_taxid, lineage}, or None if NCBI does not know it.

    """
    query = str(query).strip()
    if not query:
        return None
    if query.isdigit():
        taxid = query
    else:
        found = search(query, retmax=1)
        taxid = found[0] if found else None
    if taxid is None:
        return None
    entries = _read_entries([taxid])
    return describe(entries[0]) if entries else None


def names_for(taxids):
    """{taxid: scientific name} for several taxids in a single request."""
    taxids = sorted({str(t) for t in taxids})
    if not taxids:
        return {}
    return {str(e['TaxId']): str(e['ScientificName']) for e in _read_entries(taxids)}
