# Iterated Prisoners Dilemma

Assessing the Effectiveness of Diversity, Niceness, Forgiveness, Retaliation and Non-Enviousness as Heuristics in a Genetic Algorithm Simulation of the Iterated Prisoner’s Dilemma.
  
Group members: Eric Smith, Nik Steel, Chris Zygowski and Dwayne Alleyne.

### Installation

```sh
$ git clone https://github.com/ambi0/Iterated-Prisoners-Dilemma Iterated-Prisoners-Dilemma
$ cd Iterated-Prisoners-Dilemma
```

### To Compile and Execute

```sh
$ gcc -o genetic genetic.c
$ ./genetic < input.txt
```

### Input Structure
Each line consists of a simulation of the format:
```sh
[Population Size] [Number Genetic Iterations] [Number Compete Iterations] [Selection Rate] [Mutation Rate] [] [] [] [] [] []
```

Example input:
```sh
100 100 100 .15 .05 1 0 0 0 0 0
100 100 100 .15 .05 0 1 0 0 0 0
100 100 100 .15 .05 0 0 1 0 0 0
100 100 100 .15 .05 0 0 0 1 0 0
100 100 100 .15 .05 0 0 0 0 1 0
100 100 100 .15 .05 0 0 0 0 0 1
```
### Version
1.0.0
