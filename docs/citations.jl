# Workaround for DocumenterVitepress' DocumenterCitations extension (v0.3.5), which renders a
# `@bibliography` block as a bare ordered list. That (a) throws away the labels the citation
# style produces -- the `[1]`, `[2]`, ... of the default `:numeric` style -- and (b) numbers
# the list starting from 2 (`bullet(i) = "$(i+1). "` in DocumenterVitepress' writer.jl), so
# the numbers shown never line up with the in-text citation links.
#
# Instead, emit the style's own labels into an unordered list wrapped in a `.citation` div,
# which is what `docs/src/.vitepress/theme/overrides.css` styles.
import DocumenterVitepress as DV
using DocumenterCitations: BibliographyNode
using Documenter: MarkdownAST

function DV.render(
    io::IO, mime::MIME"text/plain", node::MarkdownAST.Node,
    bibliography::BibliographyNode, page, doc; kwargs...
)
    list = MarkdownAST.Node(MarkdownAST.List(:bullet, false))
    for item in bibliography.items
        entry = MarkdownAST.Node(MarkdownAST.Paragraph())
        if !isnothing(item.anchor_key)
            push!(
                entry.children,
                MarkdownAST.Node(MarkdownAST.HTMLInline("<a id='$(item.anchor_key)'></a>"))
            )
        end
        if !isnothing(item.label)
            append!(entry.children, collect(item.label.children))
            push!(entry.children, MarkdownAST.Node(MarkdownAST.Text(" ")))
        end
        append!(entry.children, collect(item.reference.children))
        push!(list.children, MarkdownAST.Node(MarkdownAST.Item()))
        push!(last(list.children).children, entry)
    end
    # The blank lines are required: they close the HTML block so that markdown-it processes
    # the list between the tags as Markdown rather than passing it through verbatim.
    println(io, "\n<div class=\"citation\">\n")
    DV.render(io, mime, list, list.element, page, doc; kwargs...)
    println(io, "\n</div>")
end
