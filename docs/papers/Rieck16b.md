# ‘Shall I compare thee to a network?’—Visualizing the Topological Structure of Shakespeare’s Plays

Here, we briefly document how to reproduce the results reported in the paper [‘Shall I compare thee to a network?’—Visualizing the Topological Structure of Shakespeare’s Plays](http://bastian.rieck.ru/research/Vis2016.pdf),
which was presented at the Workshop on Visualization for the Digital Humanities at IEEE VIS 2016.

While Aleph was not originally used to obtain these results, we still
want to demonstrate how to perform the analysis of the paper using
Aleph's various tools.

## Setup

First, you need to obtain all the co-occurrence networks. 

## Creating the persistence diagrams

In order to create the persistence diagrams as described in the paper,
please use the following script:

{% highlight zsh %}
#!/usr/bin/env zsh

ALEPH_DIR=$HOME/Projects/Aleph
ALEPH_BUILD_DIR=$ALEPH_DIR/build

NETWORK_ANALYSIS_BINARY=$ALEPH_BUILD_DIR/examples/network_analysis
SPLIT_OUTPUT_BINARY=$ALEPH_DIR/utilities/split_output.py

OUTPUT_DIR=/tmp

for file in $argv; do
  OUTPUT=${file:t:r}".txt"
  $NETWORK_ANALYSIS_BINARY $file 2  > $OUTPUT_DIR/$OUTPUT
  $SPLIT_OUTPUT_BINARY --prefix=d --digits=1 $OUTPUT_DIR/$OUTPUT
done
{% endhighlight %}
