def fill(array, list):
    if len(array) != len(list):
        return False
    for i in range(len(array)):
      #  if int(list[i]) >= 0:
        array[i] = int(list[i])
       # else:
       #     return False
    return True

def readConstrains(prices, costs, plan, restrictions, p):
    flag = True
    f = open("C:/Users/Markk/Downloads/Telegram Desktop/Simplex/requirements.txt")
    for i in f.readlines():
        temp = i.replace(' ', '').split(':')
        if i[0] == '#' or len(i) == 1:
            continue
        if temp[0] == 'Prices':
            flag = fill(prices, temp[1].split(','))
        if temp[0] == 'Costs':
            flag = fill(costs, temp[1].split(','))
        if temp[0] == 'Plan':
            flag = fill(plan, temp[1].split(','))
        if temp[0] == 'Restrictions':
            flag = fill(restrictions, temp[1].split(','))
        if temp[0] == 'Product1requirements':
            flag = fill(p[0], temp[1].split(','))
        if temp[0] == 'Product2requirements':
            flag = fill(p[1], temp[1].split(','))
        if temp[0] == 'Product3requirements':
            flag = fill(p[2], temp[1].split(','))
        if not flag:
            break
    return flag

def fillMatrices():
    prices = [-1, -1, -1]
    plan = [-1, -1, -1]
    restrictions = [-1, -1]
    p = [[-1, -1], [-1, -1], [-1, -1]]
    costs = [-1, -1]
    if not readConstrains(prices, costs, plan, restrictions, p):
        print("errrrrorrr")
        return None, None, None
    prices = [prices[i] - costs[0] * p[i][0] - costs[1] * p[i][1] for i in range(3)]
    A = []
    c = [-i for i in prices]
    b = restrictions + plan
    for i in range(len(b)):
        A.append([])
        A[i] = [0 for _ in range(len(p))]
    A[0] = [p[i][0] for i in range(len(p))]
    A[1] = [p[i][1] for i in range(len(p))]
    A[2] = [1, 0, 0]
    A[3] = [0, 1, 0]
    A[4] = [0, 0, 1]
    return A, b, c
