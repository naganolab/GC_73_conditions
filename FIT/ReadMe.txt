2023.09.26_create.new.genes.file.R
   input: 230728_500_genes_id_conversion.tsv
   output: genes
   Note: create the vector of gene names for FIT

2023.09.26_run.FIT.Kos.R
   input: genes
          weather.01
          attribute.all.K
          weight.all.K
          log.rpm.K
   output: Kos.result/model.K.fld.i [i, 1:400]
           Kos.result/model.K.GC.i [i, 1:400]
           Kos.result/model.K.mix.i [i, 1:400]
           Kos.result/trn.smpl.K.fld
           Kos.result/trn.smpl.K.GC
           Kos.result/trn.smpl.K.mix
   Note: FIT model parameterization
2023.09.26_run.FIT.Tak.R
   input: genes
          weather.01
          attribute.all.T
          weight.all.T
          log.rpm.T
   output: Tak.result/model.T.fld.i [i, 1:400]
           Tak.result/model.T.GC.i [i, 1:400]
           Tak.result/model.T.mix.i [i, 1:400]
           Tak.result/trn.smpl.T.fld
           Tak.result/trn.smpl.T.GC
           Tak.result/trn.smpl.T.mix
2024.09.12_process.models.K.R
   input: genes
          log.rpm.K
          weather.01
          attribute.all.K
          Kos.result/model.K.fld.i [i, 1:400]
          Kos.result/model.K.GC.i [i, 1:400]
          Kos.result/model.K.mix.i [i, 1:400]
   output: Kos.result/pred.K.fld
           Kos.result/pred.K.GC
           Kos.result/pred.K.mix
2024.09.12_process.models.T.R
   input: genes
          log.rpm.T
          weather.01
          attribute.all.T
          Tak.result/model.T.fld.i [i, 1:400]
          Tak.result/model.T.GC.i [i, 1:400]
          Tak.result/model.T.mix.i [i, 1:400]
   output: Tak.result/pred.T.fld
           Tak.result/pred.T.GC
           Tak.result/pred.T.mix
2024.09.12_summay.K.R
   input: Kos.result/pred.K.fld
          Kos.result/pred.K.GC
          Kos.result/pred.K.mix
   output: Kos.result/accuracy.heatmap.K.GC.svg
           Kos.result/accuracy.heatmap.K.fld.svg
           Kos.result/accuracy.heatmap.K.mix.svg
           Kos.result/FIT.evaluate.K.svg
           Kos.result/prob.env.choice.GC_Kos
           Kos.result/prob.env.choice.fld_Kos
           Kos.result/prob.env.choice.mix_Kos
           Kos.result/env.v.choice.band.K.svg
           Kos.result/env.v.choice.K.svg
2024.09.12_summay.T.R
   input: Tak.result/pred.T.fld
          Tak.result/pred.T.GC
          Tak.result/pred.T.mix
   output: Tak.result/accuracy.heatmap.T.GC.svg
           Tak.result/accuracy.heatmap.T.fld.svg
           Tak.result/accuracy.heatmap.T.mix.svg
           Tak.result/FIT.evaluate.T.svg
           Tak.result/prob.env.choice.GC_Tak
           Tak.result/prob.env.choice.fld_Tak
           Tak.result/prob.env.choice.mix_Tak
           Tak.result/env.v.choice.band.T.svg
           Tak.result/env.v.choice.T.svg

The "log.rpm.K", "log.rpm.T", "weight.all.K", and "weight.all.T" files used in the analysis are available upon request.

