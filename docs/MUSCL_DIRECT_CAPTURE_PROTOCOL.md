# Protocolo da captura direta MUSCL

2026-09-15. Replay de 600 m autorizado pelo prompt anexado pelo usuário.
Nenhuma mudança de física, malha, limitador ou agenda temporal.

## Estado inicial em 2026-09-20 (superado pela execução abaixo)

Este é um protocolo proposto, não um relatório de execução. Existem saídas
opcionais de fluxos nas rotinas de momento e do tracer; o observador completo,
o executor com gates e o piloto maduro ainda estão pendentes. Testes das rotinas
isoladas não demonstram neutralidade do replay. Consultar
[o documento de passagem](TORNADOGENESIS_RESEARCH_HANDOFF.md) antes de continuar.

## Balanço definido antes da captura

Em cada estágio advectivo, registrar os fluxos de momento F_alpha,i usados
pelas rotinas C-grid, i=u,v e alpha=x,y,z. As reconstruções têm m²/s².
Reconstruir separadamente delta U_alpha = -dt D_alpha F_alpha, em m/s,
com as diferenças das faces nativas, fechamento periódico e paredes verticais
exatamente como nas rotinas originais. Comparar sua soma com a mudança de
velocidade efetivamente aplicada no solver e, separadamente, no rótulo LES.
O resto de atualização em ponto flutuante será registrado, não chamado difusão.

q=C U tem s^-1. Cada contribuição C delta U tem s^-1; sua integral com
dx dy M é um incremento de circulação em m²/s. A projeção adjunta equivalente
usa (D_x^T M) delta v_c - (D_y^T M) delta u_c. Não é denominada fluxo físico
de vorticidade. A separação em direções é de fluxos de momento, não uma
decomposição cinemática em transporte, inclinação e estiramento.

Para o total: delta q = C(delta U_LES + delta U_MUSCL + delta U_outros).
Outros é a soma direta dos demais estágios observados, não uma diferença de
inventários. Para o rótulo LES: delta q_LES = C(delta U_injecao + delta U_MUSCL_LES).
As mudanças efetivas do rótulo são medidas antes/depois de cada estágio;
a reconstrução por fluxos é uma leitura independente da identidade entre pontas.

Os nove blocos, seus 17 níveis e raios 1,2; 2,4; 4,2; 6 km são os da auditoria
anterior. Para cada bloco, pré-calcular discos fixos no centro inicial, móveis
na ponta final e a média das máscaras inicial/final. O termo de máscara é
I(M_b-M_a,q_a), ou I(M_b-M_a,(q_a+q_b)/2) na convenção simétrica.
Manter as duas pontas para verificar a identidade de inventário.

Guardar por passo a integral assinada, a integral do módulo espacial e a
projeção adjunta. O módulo da integral por passo é derivável sem perda desses
escalares. Guardar também a integral do módulo do campo acumulado em cada
bloco, para distinguir cancelamento intrabloco e espacial. Os demais estágios
totais são agrupados por passo; cancelamento entre eles permanece e é declarado.

## Gates e execução

Usar snapshot 93 e os 1015 dt originais, até 3300,023494 s. O primeiro bloco
(69 passos, aproximadamente 30 s) é o piloto. Só prosseguir automaticamente
após o piloto se todos os gates passarem, incluindo o snapshot original 94.

- Observador ativo/inativo: os 12 prognósticos devem ser bit a bit iguais em
  todos os passos; comparar também com todos os snapshots originais disponíveis.
- Agenda: passos consecutivos, dt idênticos, tempos com tolerância 2e-8 s.
- Reconstrução de incrementos e fechamento: erro RMS / max(RMS dos dois lados,
  1e-12 na unidade do campo) <=1e-10; reportar também erro máximo absoluto.
- Projeções e máscara: tolerância 1e-10 relativa, com escala absoluta de
  1 m²/s para cancelamentos de integrais próximas de zero.
- Comparar o rótulo nas pontas com a proveniência antiga e sua evolução direta
  com a inferida na auditoria anterior. Separar comparações de campos e integrais.
- Falha interrompe a integração; não modificar tolerâncias para passar.

A implementação proposta deve preferir as saídas de fluxo da chamada usada
na atualização. Se recalcular pela mesma rotina e estado anterior, registrar
explicitamente que isso verifica a reconstrução do incremento, não valida
independentemente a fórmula do esquema. Os argumentos opcionais de saída
expõem fluxos sem mudar a ordem das somas de tendência. A neutralidade do
observador completo ainda deve ser demonstrada no piloto.

Os gates numéricos acima condicionam a continuação da execução. Persistência
de sinal entre raios, convenções e cadências é um resultado científico, não
um gate: uma inversão de sinal não autoriza descartar a execução ou ajustar
os critérios. Antes do piloto, fixar escalas dimensionais separadas para
velocidade (m/s), rotacional (s^-1) e circulação (m²/s), sem ajustar tolerâncias
após observar o resultado.

Estimativa anterior ao piloto: cerca de 65 MB brutos para estatísticas por
passo, mais menos de 2 MB por blocos e verificações; compressão não é assumida.
Reservar pelo menos 1 GiB livre. A referência histórica é 33 min para replay
com controle; o piloto medirá o custo com a instrumentação nova. Não salvar
campos completos novos: apenas acumuladores temporários em memória e escalares.

## Implementação e validação preliminar em 2026-09-20

O observador e o executor foram implementados em `muscl_capture.py` e
`scripts/run_muscl_direct_capture.py`. Os fluxos são recalculados pelas rotinas
nativas sobre o estado anterior imutável, e comparados com mudanças efetivas.
O primeiro teste integrado detectou cancelamento na subtração de velocidades:
erro máximo 1,78e-15 m/s, relativo 2,46e-9 para um incremento muito pequeno.
Isso foi resolvido separando explicitamente duas verificações, sem alterar
tolerâncias: dt vezes a tendência contra a soma dos fluxos; e mudança efetiva
contra `fl(U + dt*tendencia) - U`. A diferença entre a atualização efetiva e
a soma dos fluxos permanece registrada como arredondamento, não difusão.
Os três testes integrados, incluindo CPU/GPU, passaram após essa correção.
Isso libera o piloto, não demonstra ainda reprodução da tempestade madura.

## Execução concluída em 2026-09-20

O piloto de 69 passos passou e liberou automaticamente os 1015 passos.
Todos os gates passaram sem alterar tolerâncias. Os 12 prognósticos foram
neutros bit a bit em cada passo e iguais aos 18 snapshots originais disponíveis.
O fechamento por passo teve erro relativo máximo 5,58e-13; a reconstrução
dos fluxos, 1,91e-16; a comparação direta/inferida do MUSCL LES, 1,06e-15.
Tempo medido: 1285,99 s; captura HDF5: 46.422.995 bytes.

Esta entrada substitui o estado pendente no início deste documento. Fluxos
foram recalculados pela rotina nativa no estado anterior, não extraídos de
uma chamada independente da fórmula. A interpretação e todos os resultados
estão em [MUSCL_DIRECT_CAPTURE_RESULTS.md](MUSCL_DIRECT_CAPTURE_RESULTS.md).
Manifestos locais: `outputs/muscl_direct_capture_20260920/manifest.json` e
`outputs/muscl_direct_analysis_20260920/manifest.json`. Não repetir o replay
por uma instrução histórica deste protocolo.
