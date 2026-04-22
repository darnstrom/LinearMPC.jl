push!(LOAD_PATH,"../src/")
using LinearMPC
using Documenter

# ============================================================
# Tab block preprocessor
#
# Introduces three custom fenced-code block types so that Julia
# code only needs to be written *once* per example:
#
#   ```@tab
#   # julia
#   <julia code shown in the Julia tab>
#   # python
#   <python code shown in the Python tab>
#   ```
#   (display-only: no code is executed)
#
#   ```@tabexample name
#   # julia
#   <julia code shown in the Julia tab AND executed via @example>
#   # exec-only          (optional)
#   <extra lines executed but not shown in the tab, e.g. ylims!>
#   # python
#   <python code shown in the Python tab>
#   ```
#
#   ```@tabsetup name
#   # julia
#   <julia code shown in the Julia tab AND executed via @setup>
#   # python
#   <python code shown in the Python tab>
#   ```
#
# The preprocessor expands each block into:
#   • a ```@raw html``` lang-switcher showing the code, and
#   • (for @tabexample/@tabsetup) a ```@example name``` /
#     ```@setup name``` block that runs the Julia code
#     (all lines auto-hidden so only output shows).
#
# Files are modified in-place before makedocs and restored in the
# finally block below.
# ============================================================

function _escape_html(s::AbstractString)::String
    s = replace(s, "&" => "&amp;")
    s = replace(s, "<" => "&lt;")
    s = replace(s, ">" => "&gt;")
    return String(rstrip(s))
end

function _make_lang_switcher_html(julia_code::AbstractString, python_code::AbstractString)::String
    j = _escape_html(julia_code)
    p = _escape_html(python_code)
    return """<div class="lang-switcher">
<div class="lang-switcher-tabs">
<button class="lang-switcher-tab active" data-lang="julia"><img src="../../assets/julia.svg" alt="" class="lang-icon"> Julia</button>
<button class="lang-switcher-tab" data-lang="python"><img src="../../assets/python.svg" alt="" class="lang-icon"> Python</button>
</div>
<div class="lang-switcher-content active" data-lang="julia"><pre><code class="language-julia">$(j)</code></pre></div>
<div class="lang-switcher-content" data-lang="python"><pre><code class="language-python">$(p)</code></pre></div>
</div>"""
end

# Append " # hide" to every non-blank line that doesn't already have it.
function _add_hide(code::AbstractString)::String
    return join(map(split(code, "\n")) do line
        stripped = rstrip(line)
        isempty(stripped) || endswith(stripped, "# hide") ? line : line * " # hide"
    end, "\n")
end

# Parse a tab block body into (julia_display, julia_exec_only, python) sections.
function _parse_sections(body::AbstractString)
    julia_lines, exec_lines, python_lines = String[], String[], String[]
    section = :none
    for line in split(body, "\n")
        s = strip(line)
        if     s == "# julia";     section = :julia
        elseif s == "# exec-only"; section = :exec_only
        elseif s == "# python";    section = :python
        elseif section == :julia;      push!(julia_lines, line)
        elseif section == :exec_only;  push!(exec_lines, line)
        elseif section == :python;     push!(python_lines, line)
        end
    end
    return rstrip(join(julia_lines, "\n")),
           rstrip(join(exec_lines, "\n")),
           rstrip(join(python_lines, "\n"))
end

function _expand_tab_blocks(content::String)::String
    # Helper: expand all matches of a regex using eachmatch (which yields RegexMatch
    # objects with .captures), building the result string by walking through matches.
    function _apply(str::String, pat::Regex, fn)::String
        buf = IOBuffer()
        last = firstindex(str)
        for m in eachmatch(pat, str)
            print(buf, str[last:prevind(str, m.offset)])
            print(buf, fn(m))
            last = m.offset + ncodeunits(m.match)
        end
        print(buf, str[last:end])
        return String(take!(buf))
    end

    # @tabexample: display code in tab; run via @example (all code hidden, output shown)
    content = _apply(content, r"```@tabexample (\S+)\n(.*?)\n```"s, function(m)
        julia_display, julia_exec_only, python = _parse_sections(String(m.captures[2]))
        html = _make_lang_switcher_html(julia_display, python)
        exec_code = isempty(julia_exec_only) ?
            julia_display : julia_display * "\n" * julia_exec_only
        "```@raw html\n$(html)\n```\n\n```@example $(m.captures[1])\n$(_add_hide(exec_code))\n```"
    end)

    # @tabsetup: display code in tab; run via @setup (silently, no output)
    content = _apply(content, r"```@tabsetup (\S+)\n(.*?)\n```"s, function(m)
        julia_display, _, python = _parse_sections(String(m.captures[2]))
        html = _make_lang_switcher_html(julia_display, python)
        "```@raw html\n$(html)\n```\n\n```@setup $(m.captures[1])\n$(julia_display)\n```"
    end)

    # @tab: display code in tab only (no execution)
    content = _apply(content, r"```@tab\n(.*?)\n```"s, function(m)
        julia_display, _, python = _parse_sections(String(m.captures[1]))
        html = _make_lang_switcher_html(julia_display, python)
        "```@raw html\n$(html)\n```"
    end)

    return content
end

# Preprocess markdown files in-place; restore them in the finally block below.
_docs_src = joinpath(@__DIR__, "src")
_modified_files = Dict{String,String}()
for (_root, _, _files) in walkdir(_docs_src)
    for _file in filter(f -> endswith(f, ".md"), _files)
        _path = joinpath(_root, _file)
        _original = read(_path, String)
        if occursin("@tabexample", _original) || occursin("@tabsetup", _original) || occursin("```@tab\n", _original)
            _modified_files[_path] = _original
            write(_path, _expand_tab_blocks(_original))
        end
    end
end

try

makedocs(sitename = "LinearMPC.jl",
         modules  = [LinearMPC],
         format   = Documenter.HTML(assets = ["assets/lang-switcher.css", "assets/lang-switcher.js"]),
         pages=["Home" => "index.md"
                "Manual" => ["Getting started" => ["Simple example" => "manual/simple.md",
                                                   "Models" => "manual/model.md",
                                                   "Objective function" => "manual/objective.md",
                                                   "Constraints" =>  "manual/constraints.md",
                                                   "State observer" => "manual/observer.md",
                                                   #"Code generation" =>  "manual/codegen.md",
                                           ],
                              "Beyond the basics" => ["Prestabilization" => "manual/prestab.md",
                                                      "Move blocking " => "manual/moveblock.md",
                                                      "Reference Preview" => "manual/reference_preview.md",
                                                      "Disturbance Preview" => "manual/disturbance_preview.md",
                                                      "Explict MPC" => "manual/explicit.md",
                                                      "Hybrid MPC" => "manual/hybrid.md",
                                                      "Robust MPC" => "manual/robust.md",
                                                     "Game-Theoretic MPC" => "manual/game.md",
                                                     "Linear control cost" => "manual/linear_cost.md",
                                                     #"Complexity Certification" => "manual/compcert.md",
                                                     "Solver Settings" => "manual/solver.md",
                                                    ],
                             "Benchmark" => "manual/benchmark.md"
                            ]
                "Functions" => "functions.md"
               ])

deploydocs(repo="github.com/darnstrom/LinearMPC.jl", devbranch="main", push_preview=true)

finally
    for (_path, _original) in _modified_files
        write(_path, _original)
    end
end
