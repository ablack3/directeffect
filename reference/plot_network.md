# Plot the evidence network

Draws the direct-effect network — drugs as nodes, comparative estimates
as edges — from a bare network object, before any model is fitted.
Optional visual encodings: edge thickness for precision, node size for
the number of comparisons a drug appears in, node shape for the presence
of an absolute anchor, node colour for drug class, and edge labels
showing the observed effects. Repeat comparisons of the same pair draw
as separate fanned edges, and multi-component networks render with each
component laid out separately.

## Usage

``` r
plot_network(
  de,
  weight_edges = TRUE,
  size_nodes = TRUE,
  shape_anchors = TRUE,
  label_edges = FALSE,
  drug_classes = NULL
)
```

## Arguments

- de:

  A `directeffect_network` created by
  [`direct_effect_network()`](https://ablack3.github.io/directeffect/reference/direct_effect_network.md).

- weight_edges:

  Encode each comparison's precision (`1 / std_error^2`) as edge
  thickness.

- size_nodes:

  Encode the number of comparisons a drug appears in as node size.

- shape_anchors:

  Encode the presence of an absolute anchor as node shape (triangle =
  anchored, circle = not).

- label_edges:

  Label each edge with its observed effect (log scale, 2 decimals). Off
  by default.

- drug_classes:

  Optional data frame with columns `drug` and `class` to colour nodes by
  drug class. Drugs not listed are shown as unclassified.

## Value

A ggplot (ggraph) object.

## Examples

``` r
comparisons <- data.frame(
  study_id   = c("S1", "S2", "S3"),
  target     = c("A", "A", "B"),
  comparator = c("B", "C", "C"),
  estimate   = c(0.0, 0.4, 0.4),
  std_error  = c(0.05, 0.05, 0.05)
)
de <- direct_effect_network(comparisons, effect_measure = "HR")
if (requireNamespace("ggraph", quietly = TRUE)) {
  plot_network(de)
}
```
