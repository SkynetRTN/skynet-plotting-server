from random import random
from tkinter.messagebox import NO

from numpy import real
from kdtree import Star, Star_tree, haversine


# for _ in range(10):
#     print([random()*360, (random()*180 - 90), random()*90])

# data_set = [
#     [255.72209505967228, -28.6461171784997, 37.44648112979013, 0, 0],
#     [164.4228508313348, 31.538237012727265, 19.410406900825073, 0, 0],
#     [268.86874963038537, -28.68378791270294, 3.764494626341474, 0, 0],
#     [64.9597405857706, -58.64508538931031, 58.9134092413615, 0, 0],
#     [5.956041141157371, 52.705485235504966, 56.779157471954456, 0, 0],
#     [336.2900569880232, -27.36858671945955, 80.70147435153656, 0, 0],
#     [198.96250850460518, 63.02560565262908, 53.55201009387328, 0, 0],
#     [44.32967850524375, 75.3553476673082, 24.002729619112536, 0, 0],
#     [266.01769195939124, 81.36047420901798, 76.01470081724258, 0, 0],
#     [42.892934030109, 70.68336708938821, 89.71731321278779, 0, 0],
# ]

data_set = [
    [202, 44, 0, 0, 0],
    [112, 50, 0, 0, 0],
    [312, -70, 0, 0, 0],
    [36, -28, 0, 0, 0],
    [54, 77, 0, 0, 0],
]

star_list = []
for data in data_set:
    star_list.append(Star(data[0], data[1], data[2], data[3], data[4]))

dummy = Star(20, 20, 20, 1, 1)


def nn_test():
    used_list = []
    start = star_list[round(random()*len(star_list))-1]
    tree = Star_tree(start)
    used_list.append(start)
    while len(used_list) < len(star_list):
        star = star_list[round(random()*len(star_list))-1]
        if star not in used_list:
            tree.insert(star)
            used_list.append(star)
    return [tree.nn(dummy), tree, used_list]


def nn_real_test(trials, node):
    count = 0
    for _ in range(trials):
        node_count = node
        real_stars = [
            Star(random()*360, (random()*180 - 90), random()*90, 0, 0)]
        real_tree = Star_tree(real_stars[0])
        for _ in range(node_count-1):
            temp = Star(random()*360, (random()*180 - 90), random()*90, 0, 0)
            real_tree.insert(temp)
            real_stars.append(temp)
        real_target = Star(
            random()*360, (random()*180 - 90), random()*90, 0, 0)
        min_d = None
        min_star = None
        for star in real_stars:
            dist = haversine(real_target, star)
            if min_d == None or dist < min_d:
                min_star = star
                min_d = dist
        result = real_tree.nn(real_target)
        try:
            assert result == min_star
            count += 1
        except:
            # print(real_target)
            # print(real_tree)
            # print(min_star)
            # print(result)
            pass
    return (count/trials)*100


# count = 0
# if_output = False
# trial = 1000
# for _ in range(trial):
#     test = nn_test()
#     try:
#         assert test[0] == star_list[4]
#         count += 1
#     except:
#         # if not if_output:
#         #     print(test[0])
#         #     print(test[1])
#         #     list = test[2]
#         #     if_output = True
#         pass

# print(count/trial)


# tree = Star_tree(star_list[0])
# for star in star_list[1:]:
#     tree.insert(star)
# print(tree)
# print(tree.nn(Star(51, 42, 59, 0, 0)))

print(nn_real_test(1000, 250))
