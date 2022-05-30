using GeneralizedPhase
using MAT
using HTTP

fp, _ = mktemp()
HTTP.download("https://github.com/mullerlab/generalized-phase/blob/main/data/exampleData.mat?raw=true", fp)
X = matread(fp)["exampleData"]
HTTP.download("https://github.com/mullerlab/generalized-phase/blob/main/data/exampleChannel.mat?raw=true", fp)
channel = matread(fp)["x"]
