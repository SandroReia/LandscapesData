# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""

# --- Importando bibliotecas --- #
import pandas as pd 
import collections

# -------------------------- #
# ---- Definindo funcoes --- #
# -------------------------- #

# --- Definindo funcao que traduz letras para numeros 
def traduzindo_letras_numeros(D):
    # --- Criando dicionario que armazena letras 
    L_N = collections.defaultdict(lambda : -1)
    # --- Percorrendo todos os elementos do data_frame
    for e in range(len(D)):
        # --- Percorrendo todas as letras de uma sequencia 
        for a in D.aaSeq[e]:
            # --- Verificando se a letra esta no dicionario 
            if L_N[a] == -1:
                # --- Adicionando a letra 
                L_N[a] = len(L_N)
    # --- Abrindo arquivop que armazena as sequencias 
    f1 = open("elife_numero_sequencia.txt", "w")
    f2 = open("elife_numero_sequencia_numerica.txt", "w")
    # --- Percorrendo as sequencias e armazenando resultado em arquivo 
    for e in range(len(D)):
        # --- Escrevendo dados em arquivo 
        f1.write("%s\t%s\n" % (e, D.aaSeq[e]))
        # --- Escrevendo a correspondencia entre letras e numeros 
        f2.write("%s\t%s\t%s\t%s\t%s\n" % (e, L_N[D.aaSeq[e][0]], L_N[D.aaSeq[e][1]], L_N[D.aaSeq[e][2]], L_N[D.aaSeq[e][3]]))
    # --- Fechando arquivos 
    f1.close()
    f2.close()
    # --- Abrindo arquivo que salva a correspondencia letra_numero
    f1 = open("elife_letra_numero.txt", "w")
    # --- Percorrendo letras 
    for l in L_N:
        # --- Escrevendo 
        f1.write("%s\t%s\n" % (l, L_N[l]))
    # --- Fechando arquivo
    f1.close() 
    # --- Retornando nada 
    return L_N


# --- Interrompendo a execucao do programa 
#raise(SystemExit)

# --- Lendo arquivo 
raw_data_v1 = pd.read_excel("elife-16965-supp1-v4.xlsx")
raw_data_v2 = pd.read_excel("elife-16965-supp2-v4.xlsx")

# --- Criando dataframes filtrando os dados brutos
df1 = pd.DataFrame()
df1["aaSeq"] = raw_data_v1.Variants
df1["Fit"] = raw_data_v1.Fitness

df2 = pd.DataFrame()
df2["aaSeq"] = raw_data_v2.Variants
df2["Fit"] = raw_data_v2["Imputed fitness"]

# --- Concatenando dataframes
raw_data = pd.concat([df1, df2]).reset_index(drop = True)

# --- Chamando funcao que traduz as sequencias em letras para numeros
L_N = traduzindo_letras_numeros(raw_data)

# --- Imprimindo dados --- #
f1 = open("elife_seq_number.txt", "w")
for i in range(len(raw_data)):
    f1.write("%s\t%s\n" % (i, raw_data.aaSeq[i]))
f1.close()

f1 = open("elife_sequence_fitness.txt", "w")
for i in range(len(raw_data)):
    f1.write("%s\t%s\n" % (i, raw_data["Fit"][i]))
f1.close()


