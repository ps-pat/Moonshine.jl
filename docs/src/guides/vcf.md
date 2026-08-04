```@meta
CurrentModule = Moonshine
ShareDefaultModule = true
```

# Working with VCF files
```@contents
Pages = ["vcf.md"]
Depth = 2:2
```

## Importing data from VCF files
Moonshine leverages [VCFTools](https://github.com/OpenMendel/VCFTools.jl) for compatibility with [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) encoded phased diploid data. Function VCFTools.convert_ht reads phased data from a `VCF` file and is able to returns it in the form of a `BitMatrix` ( i.e. a 2-dimensional [`Base.BitArray`](@extref)). On the other hand, [`Sample`](@ref)s can be constructed from a `BitMatrix`. The following example shows how these two functions can work together.

First, we import Moonshine as well as the `VCFTools` package.
```@repl
using VCFTools
using Moonshine
```

We download some example data into a temporary file.
```@repl
vcf_file = tempname(suffix = "_vcf-sample-example.vcf.gz")
download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz", vcf_file)
```

The data is then parsed and stored into a `BitMatrix`.
```@repl
dat = convert_ht(Bool, vcf_file, trans = true, save_snp_info = true)
```
Passing `Bool` as the first argument specify that we want `convert_ht` to return a `BitMatrix`. `trans = true` "transposes" the matrix: each column contains an haplotype. Finaly, `save_snp_info = true` returns various metada along with the haplotypes. Among other information, this gives us access to the markers' positions, which are of primary importance to us. `dat` is a simple [`Core.Tuple`](@extref) with entries 1 and 4 containing the haplotype matrix and marker positions, respectively. It is now straightforward to construct a [`Sample`](@ref):
```@repl
H = Sample(dat[1], dat[4])
```

Additional genetic parameters can be specified as keyword arguments:
```@repl
H2 = Sample(dat[1], dat[4], μ = 1e-8, ρ = 1e-8, Ne = 10_000)
```
