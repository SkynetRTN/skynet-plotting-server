import math


class Star:
    def __init__(self, ra, dec, r, pmra, pmdec) -> None:
        self.ra = ra
        self.dec = dec
        self.r = r
        self.pmra = pmra
        self.pmdec = pmdec

    def __init__(self, paras) -> None:
        self.ra = float(paras[0])
        self.dec = float(paras[1])
        self.r = float(paras[2])
        self.pmra = float(paras[3])
        self.pmdec = float(paras[4])

    # return if the other star is has smaller key value compared to current star
    def compare(self, star, dim):
        if (dim % 2 == 0):
            return star.ra < self.ra
        elif (dim % 2 == 1):
            return star.dec < self.dec

    def __eq__(self, obj):
        return isinstance(obj, Star) and self.ra == obj.ra and self.dec == obj.dec and self.r == obj.r and self.pmra == obj.pmra and self.pmdec == obj.pmdec

    def __str__(self) -> str:
        return ("[%d, %d, %d]" % (self.ra, self.dec, self.r))


def haversine(star1: Star, star2: Star):
    # approx of (theta/2)^2
    dec1 = math.radians(star1.dec)
    dec2 = math.radians(star2.dec)
    ra1 = math.radians(star1.ra)
    ra2 = math.radians(star2.ra)
    return math.sin((dec1 - dec2)/2)**2 + math.cos(dec1) * math.cos(dec2) * (math.sin((ra1 - ra2)/2))**2


def partition_dist(root_star, query, level):
    dec1 = math.radians(root_star.dec)
    dec2 = math.radians(query.dec)
    ra1 = math.radians(root_star.ra)
    ra2 = math.radians(query.ra)
    if level % 2 == 0:
        delta = abs(ra1 - ra1)
        delta = delta if delta < math.pi else (
            2*math.pi - delta)
    elif level % 2 == 1:
        delta = abs(dec1-dec2)
    if level == 0:
        delta_0 = abs(0 - ra2)
        delta_0 = delta if delta < math.pi else (
            2*math.pi - delta)
        delta = delta if delta < delta_0 else delta_0
    return (delta/2)**2


class Star_tree:
    def __init__(self, star: Star):
        self.content = star
        self.left = None
        self.right = None

    def insert(self, star: Star, level=0):
        new_tree = Star_tree(star)
        current = self
        is_left = True
        if self.content.compare(star, level):
            current = self.left
        else:
            current = self.right
            is_left = False
        if current == None:
            if is_left:
                self.left = new_tree
            else:
                self.right = new_tree
        else:
            if is_left:
                self.left.insert(star, level+1)
            else:
                self.right.insert(star, level+1)

    def findmin(self, dim: int, level=0) -> Star:
        if self.is_leaf():
            return self.content
        else:
            stars = []
            if self.left != None:
                stars.append(self.left.findmin(dim, level + 1))
            if (level % 2 != dim and self.right != None):
                stars.append(self.right.findmin(dim, level + 1))
            stars.append(self.content)
            return self.min_index(stars, dim)

    def min_index(self, stars: list[Star], dim: int) -> Star:
        min_star = stars[0]
        for star in stars:
            if min_star.compare(star, dim):
                min_star = star
        return min_star

    def delete(self, star, level=0):
        if self.content == star:
            if self.right != None:
                self.content = self.right.findmin(level, level + 1)
                self.right = self.right.delete(self.content, level + 1)
            elif self.left != None:
                self.content = self.left.findmin(level, level + 1)
                self.left = self.left.delete(self.content,
                                             level + 1)
            else:
                self = None

        elif self.content.compare(star, level):
            if self.left != None:
                # look up on the left side
                self.left = self.left.delete(star, level+1)
        else:
            if self.right != None:
                # look up on the right side
                self.right = self.right.delete(star, level+1)
        return self

    def nn(self, query) -> Star:
        best = None
        best_d = None

        def nn_recur(root: Star_tree, tree: Star_tree, query: Star, level):
            nonlocal best
            nonlocal best_d
            if tree.is_leaf():
                dist = haversine(query, tree.content)
                if best_d == None:
                    best = tree.content
                    best_d = dist
                else:
                    if dist < best_d:
                        best = tree.content
                        best_d = dist
                return
            else:
                is_left = True
                if tree.content.compare(query, level):
                    if tree.left == None:
                        is_left = False
                else:
                    if tree.right != None:
                        is_left = False
                branch = tree.left if is_left else tree.right
                other_branch = tree.right if is_left else tree.left
                nn_recur(tree, branch, query, level + 1)
                dist = haversine(query, tree.content)
                if dist < best_d:
                    best_d = dist
                    best = tree.content
                if other_branch != None:
                    dist = partition_dist(tree.content, query, level)
                    if dist < best_d:
                        nn_recur(tree, other_branch, query, level + 1)
            return

        nn_recur(None, self, query, 0)
        return [best, best_d]

    def tabs(self, n) -> str:
        result = ""
        for _ in range(n):
            result += "\t"
        return result

    def __str__(self) -> str:
        return self.toString()

    def toString(self, level: int = 0) -> str:
        if self.content == None:
            return ""
        else:
            result = self.tabs(level) + str(self.content) + "\n"
            level += 1
            if self.left != None or self.right != None:
                if self.left == None:
                    result += self.tabs(level) + "[]" + "\n"
                else:
                    result += self.left.toString(level)
                if self.right == None:
                    result += self.tabs(level) + "[]" + "\n"
                else:
                    result += self.right.toString(level)
            return result

    def is_leaf(self):
        return self.left == None and self.right == None
