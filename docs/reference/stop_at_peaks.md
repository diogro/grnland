# Stop a simulation when the population reaches a peak

Checks every \`every\` generations whether at least \`fraction\` of the
population lies within \`radius\` of one of the \`peaks\`. This
reproduces the convergence rule of the original scripts with \`shape =
"box"\`, \`radius = 0.1\`, \`fraction = 0.5\`, \`every = 10\`.

## Usage

``` r
stop_at_peaks(
  peaks,
  radius = 0.1,
  fraction = 0.5,
  every = 10,
  shape = c("box", "ball")
)
```

## Arguments

- peaks:

  Matrix or data frame of peak positions (trait columns first; extra
  columns such as \`fitness\` are ignored).

- radius:

  Half-width of the box, or radius of the ball.

- fraction:

  Fraction of the population required near the peak.

- every:

  Check interval in generations.

- shape:

  \`"box"\` or \`"ball"\`.
