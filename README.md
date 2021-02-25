# Stepik exercises - Coursera Bioinformatics III - UC Sandiego

To run all the tests

```
pytest --cov=. --cov-report=html --cov-report=term-missing
```
In this course, we learn algorithm to align different sequences

## Week1

The alignment problem is transformed into the Manhattan Tourist problem, that is, finding a path in a rectangular grid graph that passes most of the tourist sights. To find the path that passes most of the tourist sights, we use dynamic programming. Imagine that the graph contains n nodes and m edges. Each edge will have a score depends on whether the edge connect to a node with a tourist sight. When we try to connect a node with a previous node, we only take the previous node with the highest score accumulated from their previous nodes.

## Week2

Apply the Manhattan Tourist problem to align two strings. Include different score for match, mismatch and insert/delete punishment.

This week there are different problems:
1. global alignment: find an alignment with highest alignment score between two strings along the whole string
2. local alignment: find aligned substring of two strings that have the highest alignment score
3. fitting alignment: find the alignment of a string that fits to another string with the highest alignment score
4. overlap alignment: find the alignment between the suffix of a string and the prefix of another string with the highest alignment score