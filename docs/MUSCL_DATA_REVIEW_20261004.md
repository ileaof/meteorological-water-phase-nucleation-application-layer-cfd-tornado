# Revisao dos dados MUSCL em 2026-10-04

## Conclusao

Os dados sustentam uma atribuicao numerica robusta: o estagio MUSCL reduz a
circulacao assinada nos discos de 2,4 e 4,2 km, a 121,659 m, nos 1015 passos
e nas tres convencoes de mascara. Nao demonstram que o MUSCL cause fisicamente
a ausencia de tornado, nem que a perda seja dissipacao artificial.

A nova leitura refina o local da perda: o anel 1,2--4,2 km e prioritario.
O nucleo, as regioes externas e os niveis superiores nao repetem universalmente
o mesmo comportamento. A interpretacao exige contabilizar a mascara.

## Verificacoes adicionais

Nesta revisao nao houve nova simulacao. Uma revisao independente por agente
conferiu tabelas, contagens e inferencias. O script
`scripts/audit_muscl_saved_capture.py` leu o HDF5 existente e verificou:

- hashes do HDF5 e dos metadados contra o manifesto original;
- continuidade exata da agenda e dos inventarios entre blocos;
- soma x+y+z+arredondamento contra incremento MUSCL por passo;
- integral direta contra adjunta e soma de passos contra soma de blocos;
- fechamento dos inventarios total/LES com movimento de mascara;
- aninhamento das mascaras, necessario para subtrair discos e formar aneis;
- desigualdades entre modulo da integral, integral do modulo e agregacao temporal.

Resultado: PASS. Diferencas absolutas maximas: soma direcional total
5,19e-13 m2/s; LES 5,68e-14; fechamento total 8,37e-11; fechamento LES
7,28e-12. A diferenca entre somar os passos e somar os blocos foi 1,40e-9
m2/s sobre o conjunto completo de termos, mascaras e alturas. Esta verificacao
escalar adicional usa rtol=1e-10 e atol=1e-8 m2/s; nao altera os gates originais
do replay, ja aprovados. Agenda e continuidade de inventarios: diferenca zero.
Os 11 testes focados de fluxos/observador tambem passaram.

Saida: `outputs/muscl_saved_audit_20261004/audit.json`, com hash do script e
da captura. Copia compartilhavel: `docs/media/storm/muscl_direct_20260920/audit_20261004.json`.

## Nucleo versus anel

Janela 2790,253490--3300,023494 s, 600 m, z=121,659 m, mascara de ponta
final. Incrementos em m2/s:

| Parcela no disco de 1,2 km | Valor |
|---|---:|
| LES local | +1.935,54 |
| MUSCL total | +2.474,12 |
| Outros operadores agrupados | -30,92 |
| Soma dos operadores | +4.378,74 |
| Movimento da mascara | -8.099,35 |
| Mudanca do inventario | -3.720,61 |

Nessa convencao, a perda observada no nucleo nao pode ser descrita como
retirada assinada pela soma dos operadores. O movimento da selecao espacial
domina o saldo negativo. Isso nao identifica a causa fisica do movimento
nem implica erro do rastreador. Na mascara simetrica, MUSCL no mesmo raio
e -4.329,62: a atribuicao por parcela depende da identidade escolhida.

No anel 1,2--4,2 km, MUSCL=-21.785,77 e e negativo nos 1015 passos.
No anel 4,2--6 km, MUSCL=+10.168,67, apesar de 232 passos negativos.
Nao ha demonstracao de que a perda de um anel seja exatamente exportada
para o outro. As parcelas sao projecoes do operador de momento.

## Dependencia vertical da mascara

Raio 4,2 km, soma MUSCL em m2/s:

| Altura | Ponta final | Simetrica | Fixa inicial |
|---|---:|---:|---:|
| 1488,29 m | +3.200,46 | -9.024,31 | -131.039,62 |
| 1642,49 m | +33.191,09 | +19.859,40 | -125.592,72 |
| 1974,38 m | +64.536,69 | +52.301,80 | -106.623,42 |

A altura de inversao observada em ponta final nao e uma altura critica fisica.
A mascara fixa mantem soma negativa em todos os 17 niveis; as moveis mudam
de sinal em alturas distintas. Elas amostram regioes/pesos distintos, com o
disco centrado no rastreador de baixo nivel, nao no eixo de cada altura.
Relacionar esse contraste ao alinhamento exige cruzar os centros por altura,
nao simplesmente atribuir a diferenca a inclinacao do vortice.

## Correcoes de interpretacao

1. Confirmacao em 1015 passos exclui uma inversao escondida pela agregacao
   temporal nos discos de referencia, sob os pesos usados. Nao testa uma
   trajetoria de mascara atualizada em cada passo nativo, que nao foi salva.
2. `F_MUSCL-F_centrado` nao isola o limitador. A diferenca inclui reconstrucao
   e upwinding relativamente a referencia escolhida. Substitui-se a proposta
   anterior de chama-la simplesmente "correcao limitada" por "correcao
   relativa a referencia centrada".
3. Os 18 snapshots sao estados de fim de passo, nao os 1015 estados anteriores
   ao MUSCL. Uma avaliacao neles mede tendencias instantaneas adicionais;
   nao reproduz a integral capturada nem justifica quadratura sem erro estimado.

## Proximo passo minimo

Sem nova integracao, avaliar nos estados congelados o fluxo nativo, uma
referencia centrada e, caso se queira distinguir reconstrucao de upwinding,
uma referencia upwind de primeira ordem. Manter faces, advectores, bordas,
rotacional e mascaras iguais; verificar a identidade de soma antes de interpretar.
Reportar tendencias em m2/s2 e incrementos com dt em m2/s separadamente.
Nenhuma referencia e assumida como fisica exata ou inserida no solver.

Priorizar o anel 1,2--4,2 km e o contraste com o nucleo. Se a contribuicao
negativa persistir em todas as referencias, enfraquece uma explicacao exclusiva
pela correcao entre elas. Se depender fortemente da referencia, identifica
sensibilidade algebrica, nao causalidade fisica ou convergencia.

A formacao anterior a 2790 s, a resolucao marginal e a separacao entre
dinamica fisica e erro numerico permanecem em aberto. Esta revisao nao
altera o diagnostico descritivo de concentracao/alinhamento e perda de circulacao.

Fontes locais: `outputs/muscl_direct_analysis_20260920/` (tabelas por termo,
anel e inventario), `outputs/muscl_direct_capture_20260920/capture.h5`,
`docs/MUSCL_DIRECT_CAPTURE_RESULTS.md`. O relatorio anterior contem o metodo,
o protocolo e a proveniencia completa do replay de setembro.
