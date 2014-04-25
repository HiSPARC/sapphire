import re
import gzip
import tables

from CIC import CIC
from store_events import store_event

DATAFILE = '../generator.h5'

mysql_escape_sequences = {r'\0': '\x00',
                          r"\'": r"'",
                          r'\"': r'"',
                          r'\b': '\b',
                          r'\n': '\n',
                          r'\r': '\r',
                          r'\t': '\t',
                          r'\Z': '\x1a',
                          r'\\': '\\',
                          r'\%': '%',
                          r'\_': '_',
                         }


def process_dump(path):
    datafile = tables.open_file(DATAFILE, 'w')

    id = 0
    buffer = gzip.open(path)
    for line in buffer:
        match = re.match("INSERT INTO `message` VALUES (.*)", line)
        if match:
            id = process_insert(datafile, match.group(1), id)
    buffer.close()

    datafile.close()


def process_insert(datafile, s, id):
    insert_pattern = re.compile(r"""(\((?:[^'\)]+|'(?:\\?.)*?')+\))""")
    value_pattern = re.compile(r"""([0-9]+|'(?:\\?.)*?')""")
    for insert in insert_pattern.finditer(s):
        values = value_pattern.findall(insert.group(1))
        event_id, type, msg = int(values[0]), int(values[2]),\
                              values[4][1:-1]

        if id != 0:
            id += 1
        else:
            id = event_id

        msg = re.sub(r'\\.', unescape_mysql_string, msg)
        if id != event_id:
            raise Exception("Regexp error: %d != %d" % (id, event_id))
        event = CIC((type, msg))
        event.parseMessage()
        event = {'header': {'eventtype_uploadcode': event.uploadCode,
                            'datetime': event.datetime,
                            'nanoseconds': event.nanoseconds,
                           },
                 'datalist': event.getEventData(),
                }
        store_event(datafile, 'kascade', 601, event)
    return id


def unescape_mysql_string(matchobj):
    seq = matchobj.group(0)
    if seq in mysql_escape_sequences:
        return mysql_escape_sequences[seq]
    else:
        return seq


if __name__ == '__main__':
    process_dump('../generator-buffer.sql.gz')
