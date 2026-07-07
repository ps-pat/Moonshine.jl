set shell := ["nu", "-c"]

doc-host-location := x"${WEBSITE_LOCATION}/content/software/documentation/moonshine"

# Deploy documentation
doc-deploy:
    julia --project=docs docs/deploy.jl {{ doc-host-location }}

# Build documentation for current commit (set `mode` to `draft` to enable draft mode)
doc-make mode="": doc-check
    julia --project=docs --eval "import Pkg; Pkg.instantiate(); Pkg.precompile()"
    time julia --project=docs -- docs/make.jl {{ mode }}

# Serve documentation locally
doc-serve:
    julia --project=docs --eval "import LiveServer; LiveServer.serve(dir = \"docs/build\")"

# Check markdown files for dumb mistakes
[working-directory: "docs"]
doc-check:
    #!/usr/bin/env nu
    print $'(ansi rb)Checking(ansi reset_bold) for common mistakes...(ansi reset)'
    let mistakes = ls **/* | get name | where $it =~ '\.md$' | each { ./check-stupid-mistakes.awk $in } | compact --empty
    if ($mistakes | is-empty) {
        print $'(ansi g)No mistake found!(ansi reset)'
    } else {
        $mistakes
    }
