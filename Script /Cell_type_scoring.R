# List genes associated with S phase and G2/M phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

m <-
  CellCycleScoring(
    m,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )
