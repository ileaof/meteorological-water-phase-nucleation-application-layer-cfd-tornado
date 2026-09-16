# Balanço discreto da circulação e do rótulo LES

Definição anterior ao cálculo, 2026-09-15. Somente pós-processamento.

## Grandezas e identidade

Seja C o rotacional vertical usado pelo tracer: média das velocidades C-grid
para centros, seguida de diferenças centradas horizontais no interior e
unilaterais de primeira ordem nas bordas externas. Não substituir C por uma
decomposição cinemática contínua. Seja q = C U_LES, em s^-1, e M uma máscara
binária de disco. A integral I(M,q) = dx dy sum(M q), por nível, tem m²/s.
Não há ponderação por densidade. Integração vertical exigiria dz e daria m³/s;
esse inventário volumétrico não é circulação de um contorno horizontal.

Em um bloco [a,b] com estados e incrementos exatamente coincidentes:

- L = C sum(delta U_LES local), em s^-1, recuperável dos incrementos nativos.
- A = q_b - q_a - L, em s^-1, evolução acumulada inferida do rótulo sob MUSCL.
- I(M_b,q_b) - I(M_a,q_a) = P + H + S, em m²/s.
- P = I(M_b,L); H = I(M_b,A); S = I(M_b-M_a,q_a).

S é a contribuição exata da troca de máscara nesta convenção de ponto final,
não um fluxo físico. A convenção simétrica usa M_bar=(M_a+M_b)/2 para P,H e
q_bar=(q_a+q_b)/2 para S. Ambas fecham exatamente; sua diferença mede a
dependência da partição com a convenção em blocos finitos.

A identificação de A com a evolução do rótulo segue do código v4: somente o
estágio LES injeta nesse rótulo, e somente o estágio advectivo o atualiza depois.
O arredondamento da subtração permanece. Isso é inferência por exclusão de
operadores, não uma medição independente dos fluxos MUSCL nem uma porcentagem
causal. Confirmar identidade da vorticidade total entre arquivo e replay nas
pontas, além da neutralidade arquivada entre replay observado e seu controle.

Para a circulação total, substituir q por C U e usar todos os incrementos
de operadores salvos. A soma dos seus rotacionais deve reproduzir q_b-q_a.
Os termos cinemáticos não entram nessa soma.

## Fronteira discreta recuperável

Para qualquer incremento nativo B=(B_u,B_v), com médias centradas B_uc,B_vc:

I(M,C B) = dx dy [sum((D_x^T M) B_vc) - sum((D_y^T M) B_uc)].

O lado direito é a soma por partes discreta. Com máscaras afastadas das
bordas externas, os pesos não nulos ficam na faixa do contorno rasterizado.
É a contribuição de borda ao incremento de circulação do operador. Sua unidade
é m²/s. Não é o fluxo advectivo de zeta através da superfície do cilindro.
É calculável independentemente pelo incremento de velocidade salvo, para LES
local e para cada operador total; deve concordar com a integral do rotacional.

Não há velocidades rotuladas nem fluxos por face/etapa nos snapshots v4.
Não se pode reconstruir independentemente essa expressão para A_LES, nem
separar seus fluxos horizontais/verticais, apenas com q_LES nas pontas.

## Amostragem e cancelamentos

Usar os estados v4 0,1,3,4,8,9,11,14,15,17 de cada caso. Em 600 m são as
dez pontas coincidentes com a sequência; somar todos os incrementos entre elas.
Os mesmos índices em 300 m fornecem nove blocos nominais comparáveis. Os tempos
exatos entre grades diferem e serão reportados, sem interpolação dos campos.
Regiões: discos de raios 1,2; 2,4; 4,2; 6 km, em todos os 17 níveis até 2 km,
centrados no rastreador de baixo nível; comparar máscara móvel e disco fixo no
centro inicial. O mesmo centro é usado em toda a coluna, sem seguir o eixo em altura.

Para cada campo de incremento F: guardar sum(M F), sum(M abs(F)); para S,
sum(abs(delta M q_a)). Separar a soma temporal assinada da soma dos módulos
das integrais por bloco. Os módulos ainda contêm cancelamentos entre passos
nativos dentro de cada bloco; não comparar suas razões como se fossem taxas
brutas invariantes à cadência. Nenhuma razão de magnitudes será chamada causal.

Gates: coincidência <=1e-8 s e mesma sequência de passos; campos finitos;
diferença de zeta entre arquivo/replay <=1e-12 s^-1; fechamento total relativo
<=1e-10; identidade da máscara e soma por partes com tolerância de arredondamento.
Arquivos existentes serão abertos somente para leitura. Novos artefatos terão
manifesto com hashes dos scripts, tabelas e seleções HDF5 efetivamente lidas.
