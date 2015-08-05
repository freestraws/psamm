
import heapq
from itertools import count


class _HeapEntry(object):
    """Heap entry"""

    def __init__(self, item):
        self._item = item
        self._valid = True

    @property
    def item(self):
        return self._item

    @property
    def valid(self):
        return self._valid

    def invalidate(self):
        self._item = None
        self._valid = False


class Heap(object):
    """Heap (min-heap)"""

    def __init__(self, it=[], key=None):
        self._counter = count()
        if key is None:
            key = lambda x: x
        self._key = key

        self._h = []
        self._entries = {}
        for x in it:
            entry = _HeapEntry(x)
            self._h.append((key(x), next(self._counter), entry))
            if x in self._entries:
                raise ValueError('Duplicate item in heap')
            self._entries[x] = entry

        heapq.heapify(self._h)

    def __len__(self):
        return len(self._entries)

    def __contains__(self, item):
        return item in self._entries

    def pop(self):
        while len(self._h) > 0:
            _, _, entry = heapq.heappop(self._h)
            if entry.valid:
                del self._entries[entry.item]
                return entry.item
        raise ValueError('Pop on empty heap')

    def push(self, item):
        if item in self._entries:
            raise ValueError('Item already in heap')
        entry = _HeapEntry(item)
        self._entries[item] = entry
        heapq.heappush(self._h, (self._key(item), next(self._counter), entry))

    def update(self, item):
        if item not in self._entries:
            raise ValueError('Item not in heap')
        self._entries[item].invalidate()

        entry = _HeapEntry(item)
        self._entries[item] = entry
        heapq.heappush(self._h, (self._key(item), next(self._counter), entry))
