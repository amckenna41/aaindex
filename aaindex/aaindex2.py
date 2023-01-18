#AAIndex2 

# class AAIndex2():

# get_ref_from_record()
# get_values_from_record

##AAINDEX 2 ###
###########
# class AAIndex2():

#     def __init__(self):
#         # Record.__init__(self)
#         self.key = None
#         self.desc = ''
#         self.ref = ''
#         self.authors = ''
#         self.title = ''
#         self.journal = ''
#         self.correlated = dict()
#         self.index = dict()
#         self.comment = ''
#         self.index = []
#         self.rows = dict()
#         self.cols = dict()
#         self._aaindex = dict()
#         self._aaindex = self._parse('data/aaindex2')

#     def extend(self, row):
#         self.index.append(row)

#     def _get(self, aai, aaj):
#         i = self.rows[aai]
#         j = self.cols[aaj]
#         return self.index[i][j]

#     def get(self, aai, aaj, d=None):
#         try:
#             return self._get(aai, aaj)
#         except:
#             pass
#         try:
#             return self._get(aaj, aai)
#         except:
#             return d

#     def __getitem__(self, aaij):
#         return self.get(aaij[0], aaij[1])

#     def median(self):
#         x = []
#         for y in self.index:
#             x.extend(filter(None, y))
#         x.sort()
#         if len(x) % 2 == 1:
#             return x[len(x) // 2]
#         return sum(x[len(x) // 2 - 1:len(x) // 2 + 1]) / 2.0

#     def _float_or_None(self, x):
#         if x == 'NA' or x == '-':
#             return None
#         return float(x)

#     def _parse(self, filename, quiet=True):
#         '''
#         Parse aaindex input file. `rec` must be `Record` for aaindex1 and
#         `MarixRecord` for aaindex2 and aaindex3.
#         '''
#         if not os.path.exists(filename):
#             if sys.version_info[0] < 3:
#                 from urllib import urlretrieve
#             else:
#                 from urllib.request import urlretrieve
#             url = 'ftp://ftp.genome.jp/pub/db/community/aaindex/' + os.path.split(filename)[1]
#             print('Downloading "%s"' % (url))
#             filename = urlretrieve(url, filename)[0]
#             print('Saved to "%s"' % (filename))
#         f = open(filename)

#         # current = rec()
#         lastkey = None

#         for line in f:
#             key = line[0:2]
#             if key[0] == ' ':
#                 key = lastkey

#             if key == '//':
#                 self._aaindex[self.key] = self
#                 self = rec()
#             elif key == 'H ':
#                 self.key = line[2:].strip()
#             elif key == 'R ':
#                 self.ref += line[2:]
#             elif key == 'D ':
#                 self.desc += line[2:]
#             elif key == 'A ':
#                 self.authors += line[2:]
#             elif key == 'T ':
#                 self.title += line[2:]
#             elif key == 'J ':
#                 self.journal += line[2:]
#             elif key == '* ':
#                 self.comment += line[2:]
#             elif key == 'C ':
#                 a = line[2:].split()
#                 for i in range(0, len(a), 2):
#                     self.correlated[a[i]] = float(a[i + 1])
#             # elif key == 'I ':
#             #     a = line[1:].split()
#             #     if a[0] != 'A/L':
#             #         current.extend(map(_float_or_None, a))
#             #     elif list(Record.aakeys) != [i[0] for i in a] + [i[-1] for i in a]:
#             #         print('Warning: wrong amino acid sequence for', current.key)
#             #     else:
#             #         try:
#             #             assert list(Record.aakeys[:10]) == [i[0] for i in a]
#             #             assert list(Record.aakeys[10:]) == [i[2] for i in a]
#             #         except:
#             #             print('Warning: wrong amino acid sequence for', current.key)
#             elif key == 'M ':
#                 a = line[2:].split()
#                 if a[0] == 'rows':
#                     if a[4] == 'rows':
#                         a.pop(4)
#                     assert a[3] == 'cols' and len(a) == 6
#                     i = 0
#                     for aa in a[2]:
#                         self.rows[aa] = i
#                         i += 1
#                     i = 0
#                     for aa in a[5]:
#                         self.cols[aa] = i
#                         i += 1
#                 else:
#                     self.extend(map(self._float_or_None, a))
#             elif not quiet:
#                 print('Warning: line starts with "%s"' % (key))

#             lastkey = key
        
#         return _aaindex

# aaindex2 = AAIndex2()

# # class Record:

# #     '''
# #     Amino acid index (AAindex) Record
# #     '''
# #     aakeys = 'ARNDCQEGHILKMFPSTWYV'

# #     def __init__(self):
# #         self.key = None
# #         self.desc = ''
# #         self.ref = ''
# #         self.authors = ''
# #         self.title = ''
# #         self.journal = ''
# #         self.correlated = dict()
# #         self.index = dict()
# #         self.comment = ''

# #     def extend(self, row):
# #         i = len(self.index)
# #         for x in row:
# #             self.index[self.aakeys[i]] = x
# #             i += 1

# #     def get(self, aai, aaj=None, d=None):
# #         assert aaj is None
# #         return self.index.get(aai, d)

# #     def __getitem__(self, aai):
# #         return self.get(aai)

# #     def median(self):
# #         x = sorted(filter(None, self.index.values()))
# #         half = len(x) // 2
# #         if len(x) % 2 == 1:
# #             return x[half]
# #         return (x[half - 1] + x[half]) / 2.0

# #     def __str__(self):
# #         desc = self.desc.replace('\n', ' ').strip()
# #         return '%s(%s: %s)' % (self.__class__.__name__, self.key, desc)


# # class MatrixRecord(Record):
# class MatrixRecord():

#     '''
#     Matrix record for mutation matrices or pair-wise contact potentials
#     '''

#     def __init__(self):
#         # Record.__init__(self)
#         self.key = None
#         self.desc = ''
#         self.ref = ''
#         self.authors = ''
#         self.title = ''
#         self.journal = ''
#         self.correlated = dict()
#         self.index = dict()
#         self.comment = ''
#         self.index = []
#         self.rows = dict()
#         self.cols = dict()

#     def extend(self, row):
#         self.index.append(row)

#     def _get(self, aai, aaj):
#         i = self.rows[aai]
#         j = self.cols[aaj]
#         return self.index[i][j]

#     def get(self, aai, aaj, d=None):
#         try:
#             return self._get(aai, aaj)
#         except:
#             pass
#         try:
#             return self._get(aaj, aai)
#         except:
#             return d

#     def __getitem__(self, aaij):
#         return self.get(aaij[0], aaij[1])

#     def median(self):
#         x = []
#         for y in self.index:
#             x.extend(filter(None, y))
#         x.sort()
#         if len(x) % 2 == 1:
#             return x[len(x) // 2]
#         return sum(x[len(x) // 2 - 1:len(x) // 2 + 1]) / 2.0

# _aaindex = dict()

# def _float_or_None(x):
#     if x == 'NA' or x == '-':
#         return None
#     return float(x)

# def _parse(filename, rec, quiet=True):
#     '''
#     Parse aaindex input file. `rec` must be `Record` for aaindex1 and
#     `MarixRecord` for aaindex2 and aaindex3.
#     '''
#     if not os.path.exists(filename):
#         if sys.version_info[0] < 3:
#             from urllib import urlretrieve
#         else:
#             from urllib.request import urlretrieve
#         url = 'ftp://ftp.genome.jp/pub/db/community/aaindex/' + os.path.split(filename)[1]
#         print('Downloading "%s"' % (url))
#         filename = urlretrieve(url, filename)[0]
#         print('Saved to "%s"' % (filename))
#     f = open(filename)

#     current = rec()
#     lastkey = None

#     for line in f:
#         key = line[0:2]
#         if key[0] == ' ':
#             key = lastkey

#         if key == '//':
#             _aaindex[current.key] = current
#             current = rec()
#         elif key == 'H ':
#             current.key = line[2:].strip()
#         elif key == 'R ':
#             current.ref += line[2:]
#         elif key == 'D ':
#             current.desc += line[2:]
#         elif key == 'A ':
#             current.authors += line[2:]
#         elif key == 'T ':
#             current.title += line[2:]
#         elif key == 'J ':
#             current.journal += line[2:]
#         elif key == '* ':
#             current.comment += line[2:]
#         elif key == 'C ':
#             a = line[2:].split()
#             for i in range(0, len(a), 2):
#                 current.correlated[a[i]] = float(a[i + 1])
#         # elif key == 'I ':
#         #     a = line[1:].split()
#         #     if a[0] != 'A/L':
#         #         current.extend(map(_float_or_None, a))
#         #     elif list(Record.aakeys) != [i[0] for i in a] + [i[-1] for i in a]:
#         #         print('Warning: wrong amino acid sequence for', current.key)
#         #     else:
#         #         try:
#         #             assert list(Record.aakeys[:10]) == [i[0] for i in a]
#         #             assert list(Record.aakeys[10:]) == [i[2] for i in a]
#         #         except:
#         #             print('Warning: wrong amino acid sequence for', current.key)
#         elif key == 'M ':
#             a = line[2:].split()
#             if a[0] == 'rows':
#                 if a[4] == 'rows':
#                     a.pop(4)
#                 assert a[3] == 'cols' and len(a) == 6
#                 i = 0
#                 for aa in a[2]:
#                     current.rows[aa] = i
#                     i += 1
#                 i = 0
#                 for aa in a[5]:
#                     current.cols[aa] = i
#                     i += 1
#             else:
#                 current.extend(map(_float_or_None, a))
#         elif not quiet:
#             print('Warning: line starts with "%s"' % (key))

#         lastkey = key

# # _parse('data/aaindex2', MatrixRecord)

# # print(_aaindex['KOLA920101'].journal)


