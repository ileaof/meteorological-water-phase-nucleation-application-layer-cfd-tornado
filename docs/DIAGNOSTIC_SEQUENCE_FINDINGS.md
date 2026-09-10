# Diagnóstico final da sequência 3D — 6 de setembro de 2026

**A cadeia falha na manutenção da concentração de vorticidade vertical baixa, depois de ela ter se formado.** Neste caso, existe vorticidade horizontal disponível, existe inclinação positiva durante o crescimento e existe uma ligação rotativa temporária até o primeiro nível resolvido. Isso não evolui para uma concentração baixa persistente com condensação até esse nível. A atribuição dos mecanismos é suportada pelo orçamento lagrangiano, mas não é exata: seu residual RMS é 18% da mudança de zeta.

Esta conclusão se refere ao **pai idealizado de 600 m**, integrado desde a inicialização até 3.900 s, com a configuração existente de superfície. Não é uma conclusão sobre todos os ninhos antigos, nem uma demonstração de resolução de um tornado. Foram adicionadas somente instrumentação passiva e análise. As equações, os parâmetros físicos e a participação da pressão na termodinâmica permaneceram inalterados.

## Sequência e validação

Foram salvos **131 estados 3D sincronizados e 6.414 passos**, com u/v/w nativos nas faces, as duas pressões, perfis do estado base, variáveis termodinâmicas, hidrometeoros, incrementos por operador e integrais cinemáticas. As tendências abrangem 0–2 km com halos; os estados cobrem a coluna completa. O arquivo tem 13.058 GB decimais, cerca de 12.16 GiB.

O orçamento discreto fecha: residual máximo de **8.28e-18 s^-2**, com maior razão RMS residual/mudança de zeta de **1.34e-13**. A divergência de massa calculada nas faces nativas tem máximo interior de **5.33e-17 kg m^-3 s^-1**. Isso comprova a contabilidade dos incrementos e a continuidade discreta nos instantes salvos; não transforma automaticamente o orçamento contínuo em exato.

Os testes incluem igualdade bit a bit com e sem observador em CPU/GPU, pressão de baixa memória, soluções analíticas de estiramento/inclinação, integração RK4, encerramento do rastreamento diante de perda ou ambiguidade e processamento de uma sequência sintética.

**64 testes selecionados passaram** na verificação final.

## Quando a concentração se forma e enfraquece

O componente ciclônico foi escolhido por uma regra predefinida perto de 1.800 s e acompanhado por sobreposição espacial, sem saltar para máximos independentes. O limiar de análise é zeta ≥0.003 s^-1, com conectividade de seis vizinhos.

| Evento | Tempo | Medida |
|---|---:|---|
| Primeira ligação do componente ao nível de 39.9 m | 2610 s | Ligação **resolvida e dependente do limiar**, não observação no solo |
| Pico do componente abaixo de 500 m | **2790 s** | zeta=0.007375 s^-1; z=491.7 m; w=4.122 m/s |
| Pico material da parcela central | **2850 s** | zeta=0.008146 s^-1; z=760.3 m |
| Queda persistente >20% do pico baixo do componente | **3030 s** | zeta=0.005627 s^-1; w no pico selecionado=0.490 m/s |
| Última ligação ao primeiro nível sob o limiar usado | 3540 s | Depois disso, o componente não alcança esse nível pelo critério |
| Associação espacial ambígua | 3840 s | O rastreamento inequívoco termina na saída anterior; nenhum ramo é escolhido |

O pico Euleriano muda de célula dentro do componente. Por isso, a queda de w no pico não é, isoladamente, uma trajetória de desaceleração do mesmo volume de ar. As trajetórias independentes resolvem essa distinção: a parcela central atinge seu máximo de zeta em 2.850 s e começa a perder a componente vertical **ainda subindo**, alcançando aproximadamente 1.09 km em 2.970 s. Só depois inicia a descida. Em 3.090 s sua zeta já é ligeiramente negativa.

## De onde vem a vorticidade e quais termos mudam de sinal

Foram lançadas 27 sementes determinísticas em torno do pico baixo, sem filtrar parcelas pelo resultado. A parcela central entra no trecho contínuo resolvido em **301 s**, perto da borda leste, em z≈177 m. Ela já tem vorticidade horizontal forte, (xi,eta)≈(-0.01442,+0.004794) s^-1, mas zeta≈-3.35e-10 s^-1. Portanto, não entra como um núcleo de rotação vertical intensa. A análise compara essa entrada ao perfil base; a origem anterior à entrada no domínio amostrado permanece censurada.

O perfil base nessa altura tem (xi,eta)≈(+0.001450,+0.01563) s^-1. A orientação da vorticidade na entrada já difere bastante da orientação inicial: não é possível atribuir todo o vetor entrante ao cisalhamento base intacto. Os termos anteriores à entrada resolvida não foram reconstruídos retroativamente.

As contribuições integradas da parcela central, separadas pelo pico **Euleriano rastreado** de 2.790 s, são:

| Contribuição para zeta | Até 2790 s | De 2790 a 3900 s |
|---|---:|---:|
| Estiramento axial | +0.003376 | +0.004584 |
| Inclinação | **+0.004430** | **-0.006527** |
| Dilatação | -0.000242 | -0.000569 |
| Rotacional da LES | -0.000052 | **-0.002298** |
| Restante do operador advectivo | -0.000676 | -0.001907 |
| Residual material | +0.000435 | -0.000203 |

Unidades: s^-1. Coriolis, arrasto, projeção e contornos estão discriminados nos CSVs; a tabela destaca as contribuições mais relevantes dessa parcela. A separação em 2.790 s não deve ser confundida com o máximo posterior de cada parcela. O JSON e o CSV de fases incluem todas as 27 sementes.

O padrão de sinal é comum às **27/27 parcelas**: inclinação integrada positiva antes do pico rastreado e negativa depois. Todas apresentam perda líquida de zeta na fase posterior. As medianas dessa fase são: inclinação **-0.006449 s^-1**, estiramento **+0.004975 s^-1**, LES **-0.001728 s^-1** e mudança observada **-0.006478 s^-1**. As medianas de termos distintos não devem ser somadas como se constituíssem a trajetória de uma única parcela.

**O estiramento não desaparece em toda a trajetória.** Ele continua contribuindo positivamente na integral posterior, mas passa a ser superado pela inclinação negativa, pela contribuição negativa da LES e pelo restante advectivo. A perda de zeta não equivale à eliminação de toda a vorticidade: ao fim, a parcela central ainda tem |omega_horizontal|≈0.0239 s^-1, embora zeta seja apenas 0.000461 s^-1.

Os orçamentos de xi e eta também foram preservados. Eles mostram contribuições opostas importantes de LES e arrasto; não seria correto interpretar uma pequena mudança líquida como ausência de geração ou de atuação desses termos. O rotacional da flutuabilidade atua diretamente nas componentes horizontais, que podem depois contribuir para zeta por inclinação.

## Pressão dinâmica e funil visível

O campo termodinâmico `p` permaneceu **identicamente zero**. A pressão da projeção `p_dyn` foi salva separadamente, sem ser inserida na saturação. Não houve condensado ql+qi acima de 1e-5 kg/kg no nível mais baixo dentro do componente rastreado.

No pico de 2.790 s, a anomalia dinâmica mediana do disco de 1.2 km em relação ao anel de 3–6 km, na altura do pico, é **-4.70 Pa**. O condensado do componente começa perto de **706 m**. A existência de pressão dinâmica não comprova um funil visível, e essa referência espacial não isola exclusivamente a contribuição centrífuga de um vórtice estreito. A figura correspondente mostra pressão, condensado e zeta separadamente.

## O que está comprovado, suportado e apenas sugestivo

**Comprovado nos dados e testes:** a instrumentação preserva os prognósticos nos testes de passividade; os incrementos discretos fecham; a pressão dinâmica permanece separada da termodinâmica; o componente ganha e perde intensidade baixa na janela observada; não há condensado no primeiro nível pelo critério informado. Também está medida a disponibilidade de vorticidade horizontal na entrada das parcelas.

**Suportado:** o crescimento da componente vertical depende de inclinação e estiramento positivos; sua manutenção falha depois, com reversão da contribuição da inclinação e contribuição negativa da LES, enquanto o ar continua sua trajetória pela tempestade. Há ligação baixa temporária, mas não concentração persistente de um núcleo estreito demonstrada por esta malha.

**Sugestivo ou não identificado:** a divisão causal precisa entre geometria, arrasto, baroclinicidade e efeito numérico como causa primária. O residual lagrangiano RMS é **18.1%**; a diferença entre o rotacional advectivo aplicado e a forma contínua tem razão RMS mediana de **41.1%**. Essa diferença inclui a forma conservativa em escoamento divergente e a discretização; não é medida pura de difusão numérica. Atribuir a falha a um único termo seria exceder a evidência.

As trajetórias são pouco sensíveis ao subpasso RK4 de 2 s versus 1 s: diferença máxima de 0.089 m. A cadência nominal de 30 s versus 60 s produz mediana de 15.7 m e máximo de 154.5 m entre pares válidos. Isso não elimina o erro do orçamento material nem resolve processos abaixo do primeiro centro ou abaixo da escala de malha.

## Arquivos para auditoria

- [Relatório numérico gerado automaticamente](../outputs/diagnostic_sequence_20260905/analysis/RELATORIO_FINAL.md)
- [Metadados finais da execução](../outputs/diagnostic_sequence_20260905/metadata.json)
- [Resumo numérico e contribuições](../outputs/diagnostic_sequence_20260905/analysis/summary.json)
- [Orçamento lagrangiano](../outputs/diagnostic_sequence_20260905/analysis/lagrangian_budget.csv)
- [Rastreamento do componente](../outputs/diagnostic_sequence_20260905/analysis/vortex_track.csv)
- [Trajetórias e evolução de zeta](../outputs/diagnostic_sequence_20260905/analysis/lagrangian_paths.png)
- [Pressão versus condensado](../outputs/diagnostic_sequence_20260905/analysis/pressure_vs_visible_cloud.png)
- [Procedimento, esquema e reprodução](DIAGNOSTIC_SEQUENCE.md)

O arquivo `source_snapshot.zip` preserva os fontes conferidos contra os hashes da execução, e `sequence.sha256.json` identifica o HDF5. A investigação termina aqui, no diagnóstico, sem propor mudanças físicas.
