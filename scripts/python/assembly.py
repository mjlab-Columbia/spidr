from collections import OrderedDict


class Assembly(object):
    _chromsizes = None
    _resolution = None

    def __init__(self, resolution):
        self._resolution = resolution
        self.init_offsets()

    def init_offsets(self):
        count = 0
        self._offsets = OrderedDict()
        for (chrom, size) in self._chromsizes.items():
            self._offsets[chrom] = count
            count += -(-size // self._resolution)

    def get_size(self, chrom):
        return self._chromsizes.get(chrom)

    def get_offset(self, chrom):
        return self._offsets.get(chrom)

    def get_position(self, n):
        for (chrom, offset) in reversed(self._offsets.items()):
            if n - offset >= 0:
                return (chrom, (n - offset) * self._resolution)
        raise ValueError

    def get_index(self, chrom, pos):
        read_bin = pos // self._resolution
        offset = self.get_offset(chrom)
        if offset is not None:
            return read_bin + offset


class Mm9(Assembly):

    def __init__(self, resolution):
        self._chromsizes = OrderedDict([
            ('chr1', 197195432),
            ('chr2', 181748087),
            ('chr3', 159599783),
            ('chr4', 155630120),
            ('chr5', 152537259),
            ('chr6', 149517037),
            ('chr7', 152524553),
            ('chr8', 131738871),
            ('chr9', 124076172),
            ('chr10', 129993255),
            ('chr11', 121843856),
            ('chr12', 121257530),
            ('chr13', 120284312),
            ('chr14', 125194864),
            ('chr15', 103494974),
            ('chr16', 98319150),
            ('chr17', 95272651),
            ('chr18', 90772031),
            ('chr19', 61342430),
            ('chrX', 166650296),
            ('chrY', 15902555),
            ('chrM', 16299)])
        super(Mm9, self).__init__(resolution)

    @classmethod
    def is_named(cls, name):
        return name == "mm9"


class Mm10(Assembly):

    def __init__(self, resolution):
        self._chromsizes = OrderedDict([
            ('chr1', 195471971),
            ('chr2', 182113224),
            ('chr3', 160039680),
            ('chr4', 156508116),
            ('chr5', 151834684),
            ('chr6', 149736546),
            ('chr7', 145441459),
            ('chr8', 129401213),
            ('chr9', 124595110),
            ('chr10', 130694993),
            ('chr11', 122082543),
            ('chr12', 120129022),
            ('chr13', 120421639),
            ('chr14', 124902244),
            ('chr15', 104043685),
            ('chr16', 98207768),
            ('chr17', 94987271),
            ('chr18', 90702639),
            ('chr19', 61431566),
            ('chrX', 171031299),
            ('chrY', 91744698),
            ('chrM', 16299)])
        super(Mm10, self).__init__(resolution)

    @classmethod
    def is_named(cls, name):
        return name == "mm10"


class Hg19(Assembly):

    def __init__(self, resolution):
        self._chromsizes = OrderedDict([
            ('chr1', 249250621),
            ('chr2', 243199373),
            ('chr3', 198022430),
            ('chr4', 191154276),
            ('chr5', 180915260),
            ('chr6', 171115067),
            ('chr7', 159138663),
            ('chr8', 146364022),
            ('chr9', 141213431),
            ('chr10', 135534747),
            ('chr11', 135006516),
            ('chr12', 133851895),
            ('chr13', 115169878),
            ('chr14', 107349540),
            ('chr15', 102531392),
            ('chr16', 90354753),
            ('chr17', 81195210),
            ('chr18', 78077248),
            ('chr19', 59128983),
            ('chr20', 63025520),
            ('chr21', 48129895),
            ('chr22', 51304566),
            ('chrX', 155270560),
            ('chrY', 59373566),
            ('chrM', 16571)])
        super(Hg19, self).__init__(resolution)

    @classmethod
    def is_named(cls, name):
        return name == "hg19"


class Hg38(Assembly):

    def __init__(self, resolution):
        self._chromsizes = OrderedDict([
            ('chr1', 248956422),
            ('chr2', 242193529),
            ('chr3', 198295559),
            ('chr4', 190214555),
            ('chr5', 181538259),
            ('chr6', 170805979),
            ('chr7', 159345973),
            ('chr8', 145138636),
            ('chr9', 138394717),
            ('chr10', 133797422),
            ('chr11', 135086622),
            ('chr12', 133275309),
            ('chr13', 114364328),
            ('chr14', 107043718),
            ('chr15', 101991189),
            ('chr16', 90338345),
            ('chr17', 83257441),
            ('chr18', 80373285),
            ('chr19', 58617616),
            ('chr20', 64444167),
            ('chr21', 46709983),
            ('chr22', 50818468),
            ('chrX', 156040895),
            ('chrY', 57227415),
            ('chrM', 16569)])
        super(Hg38, self).__init__(resolution)

    @classmethod
    def is_named(cls, name):
        return name == "hg38"


def build(name, resolution):
    for cls in Assembly.__subclasses__():
        if cls.is_named(name):
            return cls(resolution)
    raise ValueError
