# Stepik exercises - Coursera Bioinformatics III - UC Sandiego

To run all the tests

```
pytest --cov=. --cov-report=html --cov-report=term-missing
```
In this course, we learn algorithm to align different sequences

## Week1
The alignment problem is transformed into the Manhattan Tourist problem, that is, finding a path in a rectangular grid graph that passes most of the tourist sights. To find the path that passes most of the tourist sights, we use dynamic programming. Imagine that the graph contains n nodes and m edges. Each edge will have a score depends on whether the edge connect to a node with a tourist sight. When we try to connect a node with a previous node, we only take the previous node with the highest score accumulated from their previous nodes.