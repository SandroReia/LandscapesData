# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""

# --- Importando bibliotecas --- #
import pandas as pd 
import networkx as nx
import matplotlib.pyplot as plt
import collections
import numpy as np
import copy
import json
import matplotlib.patches as mpatches

from matplotlib.cm import get_cmap

#cmap = get_cmap('Reds')

#cmap = get_cmap('copper')
#cmap = get_cmap('inferno')
#cmap = get_cmap('YlOrRd')
#cmap = get_cmap('viridis')
#cmap = get_cmap('plasma')




plt.style.use('default')

# -------------------------- #
# ---- Definindo funcoes --- #
# -------------------------- #

# --- Definindo funcao que cria a rede
def criando_rede_de_configuracoes(raw_data):
    # --- Criando lista que armazena a rede 
    Rede_dict = {}    
    # --- Adicionando elementos no dicionario 
    for i in range(len(raw_data)):
        Rede_dict[raw_data.aaSeq[i]] = []        
    # --- Percorrendo dataframe e comparando dois a dois 
    for i in range(len(raw_data) - 1):
        for j in range((i + 1), len(raw_data)):            
            # --- Definindo o contador de iguais
            contador = 0                    
            # --- Percorrendo as letras do elemento i
            for k in range(len(raw_data.aaSeq[i])):
                if (raw_data.aaSeq[i][k] == raw_data.aaSeq[j][k]):
                    contador += 1                                                
            # --- Verificando se tem apenas 1 diferente:
            if contador == 8:               
                # --- Adicionando no dicionario 
                Rede_dict[raw_data.aaSeq[i]].append(raw_data.aaSeq[j])
                Rede_dict[raw_data.aaSeq[j]].append(raw_data.aaSeq[i])
                    
    # ---- Retornando 
    return Rede_dict


# --- Definindo funcao que determina os maximos locais 
def Determinando_Maximos_Locais(raw_data, rede):
    # --- Criando dicionario que armazena as fitness das configuracoes 
    dFit = {}
    # --- Adicionando configuracoes ao dicionario 
    for i in range(len(raw_data)):
        dFit[raw_data["aaSeq"][i]] = raw_data["median"][i]
    # --- Criando vetor que armazena as sequencias que compoe os maximos locais 
    mLocais = []
    # --- Percorrendo as sequencias da rede e verificando quais sao maximos locais 
    for seq in rede:
        # --- Assumindo que a sequencia atual eh maximo local
        mlocal = True
        # --- Percorrendo as sequencias vizinhas 
        for viz in rede[seq]:
            # --- Verificando se o fitness do vizinho eh maior que o fitness da configuracao atual
            if (dFit[viz] >= dFit[seq]):
                # --- Dizendo que nao eh maximo local 
                mlocal = False 
                # --- Interrompendo a execucao 
                continue
        # --- Verificando se a sequencia eh maximo local 
        if (mlocal):
            # --- Adicionando a sequencia ao vetor de maximos locais 
            mLocais.append(seq)
    # --- Retornando
    return dFit, mLocais 

# --- Definindo funcao que constroi a rede 
def Criando_Rede(dFit, Rfilt):    
    # --- Criando rede
    G = nx.DiGraph()
    # --- Adicionando no
    for key in dFit:
        # --- Adicionando nos
        G.add_node(key)        
    # --- Adicionando links 
    for key1 in Rfilt:
        for key2 in Rfilt[key1]:
            if dFit[key2] > dFit[key1]:
                G.add_edge(key1, key2)  
    # --- Retornando
    return G

# --- Definindo funcao que plota a nova rede
def plotando_rede_nova(wild, Rede_dict, DPaths):
    
    # --- Definindo a camada de cada no --- #
    
    # --- Definindo lista de nos a serem visitados
    visitar = collections.defaultdict(lambda : True)
    # --- Definindo dicionario que armazena as camadas dos nos 
    camadas = collections.defaultdict(lambda : [])
    # --- Definindo vetor de nos recem adicionados 
    recem_add = []
    # --- Definindo variavel que armazena a camada de cada no
    camada_do_no = {}
    
    # --- Definindo id da camada inicial
    camadaID = 0
    # --- Definindo nos da camada 0
    camadas[camadaID].append(wild)
    # --- Removendo wild type dos nos a serem visitados 
    visitar[wild] = False
    # --- Adicionando wild no vetor de nos recem adicionados 
    recem_add.append(wild)
    # --- Adicionando camanda do no 
    camada_do_no[wild] = 0
    
    # --- Percorrendo nos a serem visitados 
    while len(visitar) < len(Rede_dict):        
        # --- Atualizando o id da camada 
        camadaID += 1        
        # --- Criando vetor auxiliar
        adicionar = []        
        # --- Identificar primeiros vizinhos
        for no in recem_add:
            # --- Percorrendo os vizinhos 
            for viz in Rede_dict[no]:
                # --- Verificando se o no deve ser adicionado 
                if visitar[viz]:
                    # --- Adicionando no na camada
                    camadas[camadaID].append(viz)
                    # --- Informando que o no ja foi visitado
                    visitar[viz] = False
                    # --- Adicionando no nos adicionanr
                    adicionar.append(viz)                    
                    # --- Armanzenando a camada do no 
                    camada_do_no[viz] = camadaID
        # --- Atualizando nos a serem visitados
        recem_add = copy.copy(adicionar)
    
    # --- Definindo as posicoes dos nos --- #
    
    # --- Definindo raio padrao 
    raio = 10
    # --- Definindo valor de pi 
    pi = 3.14159265359
    # --- Definindo vetor de frequencias
    freqmax = max([dFit[i] for i in dFit])
    
    # --- Definindo dicionario de posicoes
    posicoes = {}
    # --- Definindo dicionario de cores 
    cores = {}
 
    # --- Definindo mapa de cores
    #cmap = get_cmap('rainbow')
    cmap = get_cmap('jet')
    
    # --- Percorrendo camadas
    for camada in range(len(camadas)):
        
        # --- Definindo a variacao angular
        dtheta = 2.0 * pi / len(camadas[camada])         
        # --- Definindo vetor de frequencias dos nos na camada
        freqs = [dFit[i] for i in camadas[camada]]
        # --- Fazendo o ranqueamento dos indices
        args = np.argsort(freqs)        
        
        # --- Percorrendo os nos da camada 
        for i in range(len(args)):            
            # --- Definindo o argumento
            arg = args[i]            
            # --- Definindo o no 
            no = camadas[camada][arg]        
            # --- Definindo posicoes x e y
            x = camada * raio * np.cos(i * dtheta)
            y = camada * raio * np.sin(i * dtheta)            
            # --- Definindo posicao do no 
            posicoes[no] = (x, y)
            cores[no] = cmap(dFit[no] / freqmax)    
            
    # --- Desenhando rede --- #
                
    # --- Definindo rede H
    H = nx.DiGraph()
    # --- Adicionando nos 
    H.add_nodes_from(Rede_dict.keys())
    
    # --- Definindo mapa de cores            
    mapa_cor = []
    for node in H.nodes():
        mapa_cor.append(cores[node])

    # --- Carregando paleta de cores
    #cmap = get_cmap('Reds')
    cmap = get_cmap('viridis')
    
    # --- Percorrendo os maximos locais 
    for local in DPaths.maximo.unique():
        # --- Fatiando o dataframe 
        fatia1 = DPaths[DPaths["maximo"] == local].sort_values(by = "freq", ascending = True).tail(1000).reset_index(drop = True)
        # --- Copiando dataframe
        fatia2 = fatia1
        # --- Normalizando
        fatia2["freq"] = 1 - fatia2["freq"] / fatia2["freq"].max()
        
        # --- Criando rede 
        H = nx.DiGraph()
        # --- Adicionando nos 
        H.add_nodes_from(G.nodes())
        
        # --- Percorrendo as configuracoes
        for c in range(len(fatia2)):
            # --- Determinando o vetor 
            v = fatia2["config"][c]
            # --- Adicionando links 
            for k in range(len(v) - 1):
                H.add_edge(v[k], v[k + 1], color = cmap(fatia2["freq"][c]), zorder = 1 - fatia2["freq"][c])
           
        # --- Definindo vetor que nos da os tamanos dos nos 
        tnos = [] 
        for no in H.nodes():
            if no == local:
                tnos.append(120)
            else:
                tnos.append(20)
            
        # --- Pegando vetor de cores
        colors = [H[u][m]['color'] for u,m in H.edges()]   
                                          
        # --- Plotando figura 
        plt.figure(figsize = (10, 10))
        # --- Desenhando os nos 
        nx.draw_networkx_nodes(H, pos = posicoes, node_size = tnos, node_color = mapa_cor, alpha = 1.0)
        nx.draw_networkx_nodes(H, pos = posicoes, nodelist = [local], node_size = 120, edgecolors = "black", node_color = cores[local], alpha = 1.0)
        nx.draw_networkx_edges(H, pos = posicoes, edge_color = colors, arrowstyle = "->", alpha = 0.5, width = 1)
        plt.axis("off")
        #plt.legend(frameon = False, handles = patchs)
        plt.tight_layout()
        plt.show()
        plt.savefig("figrede_PROB_%s.pdf" % local)
        
    # --------------------------------------------------------------------- #
    # --- Definindo funcao que plota a rede com todos os maximos locais --- #
    # --------------------------------------------------------------------- #

    # --- Percorrendo os maximos locais 
    local = DPaths.maximo.unique()

    # --- Criando rede 
    H = nx.DiGraph()
    # --- Adicionando nos 
    H.add_nodes_from(G.nodes())
      
    # --- Definindo vetor que nos da os tamanos dos nos 
    tnos = [] 
    for no in H.nodes():
        if no in local:
            tnos.append(120)
        else:
            tnos.append(20)
                                          
    # --- Plotando figura 
    plt.figure(figsize = (10, 10))
    # --- Desenhando os nos 
    nx.draw_networkx_nodes(H, pos = posicoes, node_size = tnos, node_color = mapa_cor, alpha = 1.0)
    for l in local:
        nx.draw_networkx_nodes(H, pos = posicoes, nodelist = [l], node_size = 120, edgecolors = "black", node_color = cores[l], alpha = 1.0)
    plt.axis("off")
    
    # --- Definindo as legendas 
    labels = collections.defaultdict(lambda : "None")
    # --- Definindo locais 
    locais = ["QFGWTPAME", "QFGLTALME", "QFGFSALTE", "QFGLSPLAE", "QFGLTPAQE", "QFGISALQE"]
    # --- Definindo labels 
    for i in range(len(locais)):
        labels[locais[i]] = i + 1
        
    # --- Definindo posicoes dos labels
    label_positions = copy.deepcopy(posicoes)
    # --- Definindo offset
    offset = 2.25
    for p in locais:  # raise text positions
        label_positions[p] = (posicoes[p][0] + offset, posicoes[p][1])
    p = "QFGFSALTE"
    label_positions[p] = (label_positions[p][0] - 0.25, label_positions[p][1] + 0.75)
    p = "QFGISALQE"
    label_positions[p] = (label_positions[p][0] - 0.25, label_positions[p][1] - 1)
    # --- Desenhando labels
    nx.draw_networkx_labels(H, label_positions, labels)
    
    #plt.legend(frameon = False, handles = patchs)
    plt.tight_layout()
    plt.show()
    plt.savefig("figrede_all.pdf" % local)        
    
    
    # ------------- Criando rede com links ---------------- #
    
    # --- Criando rede 
    H = nx.DiGraph()
    # --- Adicionando nos 
    H.add_nodes_from(G.nodes())
    
    # --- Plotando figura 
    plt.figure(figsize = (10, 10))
    
    # --- Desenhando os nos 
    nx.draw_networkx_nodes(H, pos = posicoes, node_size = tnos, node_color = mapa_cor, alpha = 1.0)  
    
    # --- Percorrendo os maximos locais 
    for local in DPaths.maximo.unique():
        # --- Fatiando o dataframe 
        fatia1 = DPaths[DPaths["maximo"] == local].sort_values(by = "freq", ascending = True).tail(1000).reset_index(drop = True)
        # --- Copiando dataframe
        fatia2 = fatia1
        # --- Normalizando
        fatia2["freq"] = 1 - fatia2["freq"] / fatia2["freq"].max()
                
        # --- Percorrendo as configuracoes
        for c in range(len(fatia2)):
            # --- Determinando o vetor 
            v = fatia2["config"][c]
            # --- Adicionando links 
            for k in range(len(v) - 1):
                H.add_edge(v[k], v[k + 1], color = cmap(fatia2["freq"][c]), zorder = 1 - fatia2["freq"][c])
           
        # --- Definindo vetor que nos da os tamanos dos nos 
        tnos = [] 
        for no in H.nodes():
            if no == local:
                tnos.append(120)
            else:
                tnos.append(20)
            
        # --- Pegando vetor de cores
        colors = [H[u][m]['color'] for u,m in H.edges()]   
        
        nx.draw_networkx_nodes(H, pos = posicoes, nodelist = [local], node_size = 120, edgecolors = "black", node_color = cores[local], alpha = 1.0)
        nx.draw_networkx_edges(H, pos = posicoes, edge_color = "orange", arrowstyle = "->", alpha = 0.5, width = 1)
    # --- Desenhando labels
    nx.draw_networkx_labels(H, label_positions, labels)
    plt.axis("off")
    #plt.legend(frameon = False, handles = patchs)
    plt.tight_layout()
    plt.show()
    plt.savefig("figrede_all_w_links.pdf")
        
    # --- Retornando
    return None

# --- Definindo funcao que plota a distribuicao de fitness 
def plotando_distribuicao_fitness(raw_data, title):
    # --- Calculando histograma 
    hist, edges = np.histogram(raw_data["median"], bins = 40)
    # --- Definindo centros
    centros = [((edges[i] + edges[i + 1]) / 2.0) for i in range(len(edges) -  1)]
    # --- Definindo largura das barras 
    dw = edges[1] - edges[0]
    # --- Normalizando hist 
    hist = [(hist[i]) / float(len(raw_data)) for i in range(len(hist))]
    # --- Abrindo ambiente de figura
    plt.figure()
    # --- Plotando o histograma 
    plt.bar(centros, hist, color = "orange", edgecolor = "black", width = dw, alpha = 0.6)
    # --- Plotando linha vertical 
    plt.plot([1, 1], [0, 1], linestyle = "--", color = "black", label = "wild")
    plt.xlabel("fitness")
    plt.ylabel("Frequency")
    plt.title("%s" % title)
    plt.ylim(0, 0.5)
    #plt.legend(frameon = False)
    plt.show()
    plt.tight_layout()
    plt.savefig("dist_%s.pdf" % title)
    # --- Retornando 
    return None

# ---------------------------------------------------------------------------- #

# --- Interrompendo execucao do programa
raise(SystemExit)

# --- Lendo arquivo 
raw_data = pd.read_csv("HSP90_fitness_landscape.txt", sep = "\t")

# --- Determinado configuracao wild type 
wild = raw_data["aaSeq"][0]
# --- Determinando maximo global 
mGlobal = raw_data["aaSeq"].loc[raw_data["median"].idxmax()]

# --- Lendo a rede
try: 
    Rede_dict = json.load(open("Rede_dict.txt"))
except:
    Rede_dict = criando_rede_de_configuracoes(raw_data)

# --- Criando rotina que determina os maximos locais 
dFit, mLocais = Determinando_Maximos_Locais(raw_data, Rede_dict)
            
# --- Chamando funcao que cria a rede com base no dicionario
G = Criando_Rede(dFit, Rede_dict)

# --- Definindo dicionario que armazenara os shotest paths
dPath = {}
NoPath = []
YesPath = []
# --- Calculando menores caminhos da wild type com os maximos locais 
for local in mLocais :
    # --- Calculando o shortest path entre a wild e o maximo local
    try:
        vetor = nx.shortest_path(G, source = wild, target = local)
        dPath[local] = vetor
        YesPath.append(local)
    except nx.NetworkXNoPath:
        NoPath.append(local)

try:
    AllPaths = json.load(open("AllPaths_menor+%s.txt" % 20))
except: 
    # --- Caminho a ser adicionando no maximo 
    plus_local = 20    
    # --- Definindo dicionario que armazena todos os caminhos entre duas configuracoes 
    AllPaths = collections.defaultdict(lambda : [])
    # --- Percorrendo os maximos locais 
    for local in YesPath:
        # --- Calculando todos os caminhos 
        for p in nx.all_simple_paths(G, source = wild, target = local, cutoff = len(dPath[local]) + (plus_local - 1)):
            # --- Adicionando caminho 
            AllPaths[local].append(p)            
    # --- Salvando em dicionario
    json.dump(AllPaths, open("AllPaths_menor+%s.txt" % plus_local,'w'))

# --- Criando dataframe 
Dalp = pd.DataFrame(columns = ["maximo", "len", "freq", "config"])
# --- Passando AllPaths para dataframe 
for local in AllPaths.keys():
    for i in range(len(AllPaths[local])):
        # --- Adicionando no dataframe
        Dalp.loc[len(Dalp)] = [local, int(len(AllPaths[local][i])), 1, AllPaths[local][i]]
# --- Passando para numerico 
Dalp.len = pd.to_numeric(Dalp.len)
 
# --- Chamando funcao que plota a rede nova
plotando_rede_nova(wild, Rede_dict, Dalp)

# --- Chamando funcao que plta a distribuicao
plotando_distribuicao_fitness(pd.read_csv("HSP90_fitness_landscape.txt", sep = "\t"), "HSP90")
