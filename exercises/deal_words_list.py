# 1:文字数多い単語、上位10個
def many_ten(words):

    #上から10個をリスト化
    mojisuu_max=0
    max_ten = []
    for i in words:
        if len(max_ten) < 10:
            max_ten.append(i)
        else:
            max_ten = sorted(max_ten, key=lambda x: len(x))
            if len(i) > len(max_ten[0]):
                max_ten[0] = i
        if len(i) >= mojisuu_max:
            mojisuu_max = len(i)

    #降順に並び替え
    max_ten = sorted(max_ten, key=lambda x: len(x), reverse=True)

    #同順位をまとめる
    joui=0
    rank = {"rank1":[], "rank2":[], "rank3":[], "rank4":[], "rank5":[], "rank6":[], "rank7":[], "rank8":[], "rank9":[], "rank10":[]}
    keys = ["rank1", "rank2", "rank3", "rank4", "rank5", "rank6", "rank7", "rank8", "rank9", "rank10"]

    for n in keys:
        for i in max_ten:
            if len(i) == len(max_ten[joui]):
                rank[n].append(i)
        joui += len(rank[n])
        rank[n].sort()
        if joui > 9:
            ten = n
            break

    #第10位が同率の場合、"words.txt"を再参照して追加、重複は除く
    for i in words:
        if len(i) == len(max_ten[9]) and not(i in rank[ten]):
            rank[ten].append(i)

    #順位ラベル名変更
    label2 = "rank" + str(1+len(rank["rank1"]))
    label3 = "rank" + str(1+len(rank["rank1"])+len(rank["rank2"]))
    #label4 = "rank" + str(1+len(rank["rank1"])+len(rank["rank2"])+len(rank["rank3"]))
    #label5 = "rank" + str(1+len(rank["rank1"])+len(rank["rank2"])+len(rank["rank3"])+len(rank["rank4"]))

    #空のリストを削除
    for n in keys:
        if not(len(rank[n])):
            rank.pop(n)

    #結果表示
    print("文字数が多い単語ランキングは")
    print("rank1 :",len(rank["rank1"][0]),"characters ; ",rank["rank1"])
    print(label2,":",len(rank["rank2"][0]),"characters ; ",rank["rank2"])
    print(label3,":",len(rank["rank3"][0]),"characters ; ",rank["rank3"])
    #print(label4,":",len(rank["rank4"][0]),"characters ; ",rank["rank4"])
    #print(label5,":",len(rank["rank5"][0]),"characters ; ",rank["rank5"])



# 2:文字数の最頻値
def saihinchi(words):
    #各文字数の単語数を得る
    mojisuu_max=21
    mojisuu = [0]
    for n in range(1,mojisuu_max+1):
        cnt = 0
        for i in words:
            if len(i) == n:
                cnt += 1
        mojisuu.append(cnt)
    #[print(n,"characters : ",mojisuu[n], "words") for n in range(1,mojisuu_max+1)]

    #最多出現文字を決定
    char_max= []
    for n in range(1,mojisuu_max+1):
        if mojisuu[n] == max(mojisuu):
            char_max.append(n)
    [print("文字数の最頻値は",char_max[i],"characters ; ",mojisuu[char_max[i]],"words") for i in range(len(char_max))]



# 3:各文字の出現回数
def alphabet_cnt(words):
    from collections import defaultdict
    counter = defaultdict(lambda: 0)

    for n in words:
        for c in n:
            counter[c] += 1
    print("各文字の出現回数は")
    for k,v in sorted(counter.items()):
        print(str(k) + ": " + str(v),"characters")



# 4:同じ文字を含まない単語で、文字数最大のもの
def most_no_duplication(words):
    max = 0
    words1=[]
    words2=[]

    #100語までの暫定の最大文字数でフィルターをかける
    for n in words[0:100]:
        if len(n) > max:
            max = len(n)
    for n in words:
        if len(n) >= max:
            words1.append(n)

    #文字が重複するものを除く        
    for n in words1:
        cnt=0
        for c in n:
            if n.count(c) == 1:
                cnt +=1
        if cnt == len(n):
            words2.append(n)

    #文字数最大のものを決定
    mojisuu_max=0
    max_one_sub = []
    max_one = []

    for n in words2:
        if len(n) >= mojisuu_max:
            mojisuu_max = len(n)
            max_one_sub.append(n)
    for n in max_one_sub:
        if len(n) == mojisuu_max:
            max_one.append(n)
    
    [print("文字の重複がない文字数最大のものは", max_one[n], ";", len(max_one[n]), "characters") for n in range(len(max_one))]




if __name__ == "__main__":
    #"words.txt"から単語を読み、リスト化
    words = [line.strip() for line in open("words.txt").readlines()]

    many_ten(words)
    print("""
""")
    saihinchi(words)
    print("""
""")
    alphabet_cnt(words)
    print("""
""")
    most_no_duplication(words)
