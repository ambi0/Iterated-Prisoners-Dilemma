/* Iterated Prisoner's Dilemma Genetic Algorithm
 *
 * Assessing the Effectiveness of Diversity, Niceness, Forgiveness,
 * Retaliation and Non-Enviousness as Heuristics in a Genetic Algorithm
 * Simulation of the Iterated Prisoner’s Dilemma.
 *
 * Group members: Eric Smith, Nik Steel, Chris Zygowski and Dwayne Alleyne
 *
 * Copyright (c) 2015
 */

#include <cassert>
#include <bitset>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <sys/stat.h>


int min(int i, int j)
{
   if (i < j)
      return i;
   return j;
}

int max(int i, int j)
{
   if (j < i)
      return i;
   return j;
}

static constexpr size_t pairCombinations(size_t n) { return n * (n-1) / 2; }

//return true if and only if the file exists, false else
bool fileExists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

using Stat_array = std::array<double, 6>;

//functions for random number generation
class Rand
{
   public:
      static unsigned long rand64() 
      { 
         static std::random_device rd;
         static std::mt19937 gen(rd());
         static std::uniform_int_distribution<> dis(0, 63);
         return dis(gen);
      }
      
		static unsigned long rand63plus1() 
      { 
         static std::random_device rd;
         static std::mt19937 gen(rd());
         static std::uniform_int_distribution<> dis(1, 63);
         return dis(gen);
      }
		
      static unsigned long randBit() 
      { 
         static std::random_device rd;
         static std::mt19937 gen(rd());
         static std::uniform_int_distribution<> dis(0, 1);
         return dis(gen);
      }
           
      static unsigned long randBinN(unsigned int n, double p)
      {
         static std::random_device rd;
         static std::mt19937 gen(rd());
         std::binomial_distribution<> dis(n, p);
         return dis(gen);
      }
      
      static unsigned long randn(unsigned int n)
      {
         static std::random_device rd;
         static std::mt19937 gen(rd());
         std::uniform_int_distribution<> dis(0,n-1);
         return dis(gen);
      }
      
      static std::bitset<64> randBitset64()
      {
         std::bitset<64> bits = 0;
         for (auto i = 0; i < 64; ++i)
            bits[i] = randBit();
         return std::move(bits);
      }
      
      static std::bitset<6> randBitset6()
      {
         std::bitset<6> bits = 0;
         for (auto i = 0; i < 6; ++i)
            bits[i] = randBit();
         return std::move(bits);
      }
};

class Prisoner
{
   public:
      enum Choice { Cooperate=0, Defect=1 };
      
      class History
      {
         public:
            History() : hist(Rand::randBitset6()) {};
            History(std::string s) : hist(s) {};
            History(unsigned long l) : hist(l) {};
            
            void remember(Choice A, Choice B) {
               hist <<= 1;
               hist |= A;
               hist <<= 1;
               hist |= B;
            }
            unsigned long get() const {
               return hist.to_ulong();
            }
            
            static History invert(History const& H)
            {
               static std::bitset<6> even("101010"), odd("010101");
               return std::move(History((((H.hist&even)>>1)|((H.hist&odd)<<1)).to_ulong()));
            }
            
            //returns true if the opponent's bits are all 1's
            static bool paranoid(History const& H)
            {
               static std::bitset<6> opp("010101");
               return (H.hist & opp) == opp;
            }
            
            //returns true if any of the opponent's bits are 1's
            static bool uneasy(History const& H)
            {
               static std::bitset<6> opp("010101"), zero(0);
               return (H.hist & opp) != zero;
            }
            
            //returns true if any pair of bits is "01"
            static bool burned(History const& H)
            {
               return (   ((H.hist[0]==1) && (H.hist[1]==0))
                        ||((H.hist[2]==1) && (H.hist[3]==0))
                        ||((H.hist[4]==1) && (H.hist[5]==0)) );
            }
            
            static bool bothBurned(History const &H)
            {
               int i = 0, j = 0;
               if ((H.hist[0]==1) && (H.hist[1]==0)) ++i;
               if ((H.hist[0]==0) && (H.hist[1]==1)) ++j;
               if ((H.hist[2]==1) && (H.hist[3]==0)) ++i;
               if ((H.hist[2]==0) && (H.hist[3]==1)) ++j;
               if ((H.hist[4]==1) && (H.hist[5]==0)) ++i;
               if ((H.hist[4]==0) && (H.hist[5]==1)) ++j;
               return i == j;
            }
            
            friend std::ostream& operator<<(std::ostream& os, Prisoner::History const& h);
         private:
            std::bitset<6> hist;
      };
      class Strategy
      {
         public:
            Strategy() : strat(Rand::randBitset64()) {}
            Strategy(std::string s) : strat(s) {};
            
            static Strategy crossover(Strategy A, Strategy B)
            {
               unsigned long cpoint = Rand::rand63plus1();
               Strategy next;
               next.strat = (A.strat << cpoint) | (B.strat >> (64 - cpoint));
               return std::move(next);
            } 
            
            static unsigned long countUniqueBits(Strategy const& a, Strategy const& b)
            {
               std::bitset<64> s(a.strat ^ b.strat);
               unsigned long sum = 0;
               for (int i = 0; i < 64; ++i)
                  sum += s[i];
               return sum;
            }
            
            Choice choose(History const& H) const
            {
               return ((strat[H.get()]) ? Choice::Defect : Choice::Cooperate);
            }

            void mutate()
            {
               unsigned long mpoint = Rand::rand64();
               strat[mpoint].flip();
            }
            
            int forgiveness(History const & h) const
            {
               if ( Prisoner::History::uneasy(h) 
                    && !Prisoner::History::burned(h) 
                    && (choose(h) == Choice::Cooperate) )
                  return 1;
               if ( Prisoner::History::uneasy(h) 
                    && !Prisoner::History::burned(h) 
                    && (choose(h) == Choice::Defect) )  
                  return -1;
               //otherwise
               return 0;
            }
            
            double calc_forgiveness()
            {
               double sum;
               for (unsigned long i = 0; i < 64; ++i)
                  sum += forgiveness(History(i));  
               assert((sum >= -19) && (sum <= 19)); 
               return (sum + 19)/(19*2);
            }
            
            int niceness(History const & h) const
            {
               if ( !Prisoner::History::uneasy(h)  
                    && (choose(h) == Choice::Cooperate) )
                  return 1;
               if ( !Prisoner::History::uneasy(h) 
                    && (choose(h) == Choice::Defect) )  
                  return -1;
               //otherwise
               return 0;
            }
            
            double calc_niceness()
            {
               double sum;
               for (unsigned long i = 0; i < 64; ++i)
                  sum += niceness(History(i));   
               assert((sum >= -8) && (sum <= 8));
               return (sum + 8)/(8*2);
            }
            
            int retaliation(History const & h) const
            {
               if ( ( Prisoner::History::burned(h) 
                      || Prisoner::History::paranoid(h) )
                    && (choose(h) == Choice::Defect)    )
                  return 1;
               if ( ( Prisoner::History::burned(h) 
                      || Prisoner::History::paranoid(h) )
                    && (choose(h) == Choice::Cooperate) ) 
                  return -1;
               //otherwise
               return 0;
            }
            
            double calc_retaliation()
            {
               double sum;
               for (unsigned long i = 0; i < 64; ++i)
                  sum += retaliation(History(i)); 
               assert((sum >= -38) && (sum <= 38));
               return (sum + 38)/(38*2);
            }
            
            int nonEnvious(History const & h) const
            {
               if ( Prisoner::History::bothBurned(h) 
                    && (choose(h) == Choice::Cooperate) )
                  return 1;
               if ( Prisoner::History::bothBurned(h) 
                    && (choose(h) == Choice::Defect) )  
                  return -1;
               //otherwise
               return 0;
            }
            
            double calc_nonenvy()
            {
               double sum;
               for (unsigned long i = 0; i < 64; ++i)
                  sum += nonEnvious(History(i)); 
               assert((sum >= -20) && (sum <= 20));
               return (sum + 20)/(20*2);
            }
            
            friend std::ostream& operator<<(std::ostream& os, Strategy const& s);
         private:
            std::bitset<64> strat;
      };
      Prisoner() : strat(), fitness(7.75), score_ave(7.75), diversity_ave(32) {}
      Prisoner(Strategy _strat) : strat(_strat), fitness(7.75), score_ave(7.75), diversity_ave(32) {}
      
      Strategy strat;
      double fitness;
      double score_ave;
      double diversity_ave;
      Stat_array stats;
      
      //fitness is a real number between 0.0 and 1.0, for which higher is better
      void calc_fitness(Stat_array const& weights, double const& max_score, 
         double const& max_diversity, double const& min_score, double const& min_diversity)
      {
         //assertions to test the code for bugs
         assert(weights.size() == stats.size());
         
         //for score, a low value implies a high fitness
         if ((int(max_score*100000) != 0) && (max_score != min_score))
            stats[0] = 1.0 - ((score_ave - min_score) / (max_score - min_score)); 
         else
            stats[0] = 1.0;
         //for all others a high value implies high fitness
         if ((int(max_diversity*100000) != 0) && (min_diversity != max_diversity))
            stats[1] = (diversity_ave - min_diversity) / (max_diversity - min_diversity);
         else
            stats[1] = 1.0;
         stats[2] = strat.calc_niceness();
         stats[3] = strat.calc_forgiveness();
         stats[4] = strat.calc_retaliation();
         stats[5] = strat.calc_nonenvy();
         
         fitness = 0.0;
         
         //sum the weighted average of each statistic
         for (size_t i = 0; i < weights.size(); ++i)
            fitness += stats[i]*weights[i];
         
         //preferably, the weights should sum to 1. and fitness should be 0 to 1
         //assert((fitness >= 0.0) && (fitness <= 1.0));
      }
      
      void calc_aves(double const pop_size, double const num_games)
      {
         score_ave = double(score_sum) / (pairCombinations(pop_size) * num_games);
         diversity_ave = double(diversity_sum) / pairCombinations(pop_size);
      }
      
      void clear_sums()
      {
         score_sum = 0;
         diversity_sum = 0;
      }
                  
      //count the number of bits that differ between two strategies
      static void diversity(Prisoner & a, Prisoner& b)
      {
         unsigned long sum = Strategy::countUniqueBits(a.strat, b.strat);
         a.diversity_sum += sum;
         b.diversity_sum += sum;
      }
      
      static unsigned long score(Choice A, Choice B)
      {
         if (A == Choice::Cooperate && B == Choice::Cooperate)
            return 1;
         if (A == Choice::Cooperate && B == Choice::Defect)
            return 20;
         if (A == Choice::Defect && B == Choice::Cooperate)
            return 0;
         //otherwise, (A == Choice::Defect && B == Choice::Defect)
            return 10;
      }
      
      enum Result_t
      {
         allDefect, allCooperate, mixed
      };         
      
      static Result_t compete(Prisoner & A, Prisoner & B, History & H)
      {
         Choice a = A.strat.choose(H);
         Choice b = B.strat.choose(H);
         A.score_sum += score(a,b);
         B.score_sum += score(b,a);
         H.remember(a,b);
         
         if (a==Cooperate && b==Cooperate)
            return allCooperate;
         if (a==Defect && b==Defect)
            return allDefect;
         //otherwise
         return mixed;
      }
      
      friend std::ostream& operator<<(std::ostream& os, Prisoner const& p);
   private:
      unsigned long score_sum;
      unsigned long diversity_sum;
};

class Population
{
   public:
      std::vector<Prisoner> agents;
      
      //n specifies the number of agents in the vector to default construct
      Population(size_t n) : agents(n) {}
      Population(Population& p) = default;
      Population(Population&& p) = default;
      
      //sort the population by fitness in ascending order
      void sort()
      {
         std::sort(std::begin(agents),std::end(agents),
            [](Prisoner const& A, Prisoner const& B)->bool {return B.fitness < A.fitness;});
      }
      
      void get_max(double& max_sco, double& max_div)
      {
         assert(agents.size()!=0);
         
         auto first = std::begin(agents); 
         auto last = std::end(agents);

         max_sco = first->score_ave;
         max_div = first->diversity_ave;
         
         while (++first!=last)
         {
            if (max_sco < first->score_ave)
               max_sco = first->score_ave;
            if (max_div < first->diversity_ave)
               max_div = first->diversity_ave;
         } 
      }
      
      void get_min(double& min_sco, double& min_div)
      {
         assert(agents.size()!=0);
         
         auto first = std::begin(agents); 
         auto last = std::end(agents);

         min_sco = first->score_ave;
         min_div = first->diversity_ave;
         
         while (++first!=last)
         {
            if (min_sco > first->score_ave)
               min_sco = first->score_ave;
            if (min_div > first->diversity_ave)
               min_div = first->diversity_ave;
         } 
      }      
      Prisoner* get_random_agent()
      {
         if (size() < 2) return &*std::begin(agents);
         //otherwise
         return &*(std::begin(agents)+Rand::randn(size()));
      }
      
      double fitness_mean() const
      {
         double sum = 0.0;
         for (auto a : agents)
            sum += a.fitness;
         return sum / double(size());
      }
      
      double fitness_sd(double mean) const
      {
         double sum = 0.0;
         for(auto a : agents)
            sum += (a.fitness - mean) * (a.fitness - mean);
         return std::sqrt(sum / double(size()));
      }
      
      double score_mean() const
      {
         double sum = 0.0;
         for (auto a : agents)
            sum += a.score_ave;
         return sum / double(size());
      }
      
      double score_sd(double mean) const
      {
         double sum = 0.0;
         for(auto a : agents)
            sum += (a.score_ave - mean) * (a.score_ave - mean);
         return std::sqrt(sum / double(size()));
      }
      
      double diversity_mean() const
      {
         double sum = 0.0;
         for (auto a : agents)
            sum += a.diversity_ave;
         return sum / double(size());
      }
         
      double diversity_sd(double mean) const
      {
         double sum = 0.0;
         for(auto a : agents)
            sum += (a.diversity_ave - mean) * (a.diversity_ave - mean);
         return std::sqrt(sum / double(size()));
      }
      
      double stat_mean(size_t i) const
      {
         double sum = 0.0;
         for (auto a : agents)
            sum += a.stats[i];
         return sum / double(size());
      }
      
      size_t size() const { return agents.size(); }
      
      friend std::ostream& operator<<(std::ostream& os, Population const& p);
};

class Simulation
{
   public:
      struct Parameters
      {
         unsigned int pop_size;
         size_t num_genetic_iter;
         size_t num_compete_iter;
         double selection_rate;
         double mutation_rate;
         Stat_array weights;
         
         friend std::ostream& operator<<(std::ostream& os, Simulation::Parameters const& p);
      };

      //open output files, ensuring unique names
      Simulation(Parameters _params) : params(_params)
      {
         make_csv("demographics", demographics_csv);
         print_demographics_headings();
         make_csv("finalpop", finalpop_csv);
         print_finalpop_headings();
         make_csv("allpop", allpop_csv);
      }
      
      //make a file and output the simulation parameters to it
      void make_csv(std::string s, std::ofstream& os)
      {
         static unsigned long i = 0;
         while (fileExists(s + std::to_string(i) + ".csv"))
            ++i;
         os.open(s + std::to_string(i) + ".csv",
            std::ofstream::out | std::ofstream::trunc);
         os << params;
      }
      
      //close output files
      ~Simulation()
      {
         demographics_csv.close();
         finalpop_csv.close();
         allpop_csv.close();
      }
      
      //column headings
      void print_demographics_headings()
      {
         demographics_csv
            << std::endl
            << "Generation, Selected, "
            << "Fitness Mean, Fitness SD, "
            << "Score Mean, Score SD, "
            << "Diversity Mean, Diversity SD, "
            << "Relative Score Mean, Relative Diversity Mean, "
            << "Niceness Mean, Forgiveness Mean, "
            << "Retaliation Mean, Non-Envy Mean, "
            << "All-Coop Prop, All-Defect Prop, Mixed Prop"
            << std::endl;
      }
      
      //for file logging of demographics
      void log_demographics(size_t i, Population const& p)
      {
         //iterations number
         demographics_csv
            << i << ", " << num_selected << ", ";
            
         //a temporary value
         double tmp; 
         
         //output fitness average and standard deviation
         tmp = p.fitness_mean();
         demographics_csv
            << tmp << ", " << p.fitness_sd(tmp) << ", ";
            
         //output score average and standard deviation
         tmp = p.score_mean();
         demographics_csv
            << tmp << ", " << p.score_sd(tmp) << ", ";
            
         //output the diversity average and standard deviation
         tmp = p.diversity_mean();
         demographics_csv
            << tmp << ", " << p.diversity_sd(tmp);
            
         //output the six fitness parameters
         for (size_t a = 0; a < params.weights.size(); ++a)
            demographics_csv << ", " << p.stat_mean(a);
         
         //output the proportion of each move type
         demographics_csv << ", "
            << coop_prop << ", "
            << defect_prop << ", "
            << mixed_prop << std::endl;
      }
      
      void print_finalpop_headings()
      {
         finalpop_csv
            << std::endl
            << "Rank, Strategy, Fitness, Average Score, Average Diversity, "
            << "Relative Average Score, Relative Average Diversity, "
            << "Niceness, Forgiveness, Retaliation, Non-Envy"
            << std::endl;
      }
      
      //export the final population to a file in csv format
      //assumed that this function will only be run once
      void log_finalpop(Population const& p)
      {
         size_t rank = 0;
         for (auto a : p.agents)
         {
            finalpop_csv
               << ++rank << ", "
               << a.strat << ", "
               << a.fitness << ", "
               << a.score_ave << ", "
               << a.diversity_ave;
            for (auto b : a.stats)
               finalpop_csv << ", " << b;
            finalpop_csv << std::endl;
         }
      }
      
      //run the iterated prisoner's dilemma game on each unique pair
      //in the population and then calculate the fitness of each member
      void evaluate(Population &p)
      {
         //clear each player's score sum and diversity sum
         for (auto& a : p.agents)
            a.clear_sums();
         
         //clear the tally of choice types
         clear_choice_tally();
         
         //each unique pair combination of prisoners compete and get scored
         for(auto A = std::begin(p.agents); A != std::end(p.agents); ++A) {
            for(auto B = next(A); B != std::end(p.agents); ++B) {
               Prisoner::History H;
               for (size_t i = 0; i < params.num_compete_iter; ++i)
                  tally_choice(Prisoner::compete(*A,*B,H));
               Prisoner::diversity(*A,*B);
            }
         }
         
         //calculate the proportion of each choice
         calc_choice_prop();
         
         //calculate each player's average score and diversity
         for (auto& a : p.agents)
            a.calc_aves(params.pop_size, params.num_compete_iter);
         
         //determine the max score and max diversity
         double max_sco, max_div, min_sco, min_div;
         p.get_max(max_sco, max_div);
         p.get_min(min_sco, min_div);
         
         //calculate each player's fitness
         for (auto& a : p.agents)
            a.calc_fitness(params.weights, max_sco, max_div, min_sco, min_div);
         
         //sort by fitness so that other GA operations may assume a sorted population
         p.sort();
      }
      
      //clear the count of each choice type
      void clear_choice_tally()
      {
         coop_sum = 0;
         defect_sum = 0;
         mixed_sum = 0;
      }
      
      //count which choice type was made
      void tally_choice(Prisoner::Result_t r)
      {
         switch (r) 
         {
            case Prisoner::Result_t::allCooperate: ++coop_sum; break;
            case Prisoner::Result_t::allDefect: ++defect_sum; break;
            case Prisoner::Result_t::mixed: ++mixed_sum; break;
         }
      }
      
      //calculate the proportion of choices
      void calc_choice_prop()
      {
         double total = coop_sum + defect_sum + mixed_sum;
         coop_prop = coop_sum / total;
         defect_prop = defect_sum / total;
         mixed_prop = mixed_sum / total;
      }
      
      //selects a specified number of prisoners with the best fitness
      //and constructs a new population composed of them      
      //assumes the source population is sorted by fitness
      Population select_best(Population const& a)
      {
         //construct a new (empty) population
         Population b(0);
         
         //construct new prisoners using the strategies of those from the source
         for (size_t i = 0; i < num_selected; ++i)
            b.agents.emplace_back(a.agents[i].strat);
         
         //return the new population
         return std::move(b);
      }
      
      //applies crossover to all combinations of population members
      Population recombine(Population const& p)
      {
         Population parents(select_best(p));
         Population children(0);
         //for each unique pair combination of the current population
         //create a child using the crossover of the parents
         for(auto A = std::begin(parents.agents); A != std::end(parents.agents); ++A) {
            for(auto B = next(A); B != std::end(parents.agents); ++B) {
               if (children.size() < params.pop_size)
                  children.agents.emplace_back(
                     Prisoner::Strategy::crossover(A->strat,B->strat));
               else 
                  return std::move(children);
            }
         }
         return std::move(children);
      }
      
      //select a proportion of the population to apply mutation to
      void mutate(Population& p)
      {
         for (size_t i = 0; i < floor(p.size()*params.mutation_rate); ++i)
            p.get_random_agent()->strat.mutate();
      }
      
      //sorts both population, and for each high scoring element of the source population p
      //a low-scoring element from the destination is deleted and replaced
      //assumes both populations are sorted by fitness
      void replace_worst(Population &dest, Population& source)
      {
         Population next(select_best(source));
         
         //sanity check that the destination population is larger than the source
         assert(dest.size() >= next.size());
         
         // for each member of source, delete the worst from destination 
         // population
         assert(num_selected == next.size());
         for (size_t i=0; i < num_selected; ++i)
            dest.agents.pop_back();
         
         //replace with the best from source population
         for (auto& best : next.agents)
            dest.agents.emplace_back(best);
      }
      
      //the genetic algorithm
      void genetic()
      {
         Population initial(params.pop_size);
            
         for (size_t i=0; i < params.num_genetic_iter; ++i)
         {
            std::cout << "----------------------------------" << std::endl;
            std::cout << "------- Iteration " << i << " --------------" << std::endl;
            
            evaluate(initial);
         
            //determine how many to select for next pop
            num_selected = Rand::randBinN(initial.size(),params.selection_rate);
         
            //log population demographics to a csv file
            log_demographics(i, initial);
            
            std::cout << "------- Initial Population -------" << std::endl;
            std::cout << initial;
            allpop_csv << initial << std::endl;
            
            Population next = recombine(initial);
            mutate(next);
            evaluate(next);
            
            //determine how many to select for combined pop
            num_selected = Rand::randBinN(next.size(),params.selection_rate);
            
            std::cout << "------- Next Population ----------" << std::endl;
            std::cout << next;
            allpop_csv << next << std::endl;
            
            replace_worst(initial, next);
         }
         
         std::cout << std::endl << std::endl;
         
         //log the final population to a file
         log_finalpop(initial);
      }
      
   private:
      std::ofstream demographics_csv, finalpop_csv, allpop_csv;
      unsigned long defect_sum, coop_sum, mixed_sum;
      size_t num_selected;
      double defect_prop, coop_prop, mixed_prop;
      const Parameters params;
};

//I/O functions

std::ostream& operator<<(std::ostream& os, Prisoner::History const& h)
{
   os << h.hist;
   return os;
}

std::ostream& operator<<(std::ostream& os, Prisoner::Choice const& c)
{
   os << ((c == Prisoner::Choice::Cooperate) ? "cooperate" : "defect");
   return os;
}

std::ostream& operator<<(std::ostream& os, Prisoner::Strategy const& s)
{
   os << s.strat;
   return os;
}

std::ostream& operator<<(std::ostream& os, Prisoner const& p)
{
   os << p.strat << ", " << p.fitness;
   
   return os;
}

std::ostream& operator<<(std::ostream& os, Population const& p)
{
   for (auto a : p.agents)
      os << a << std::endl;
   return os;
}

std::ostream& operator<<(std::ostream& os, Simulation::Parameters const& p)
{
   os << "Population Size, " << p.pop_size << std::endl
      << "GA Iterations, " << p.num_genetic_iter << std::endl
      << "IPD Iterations, " << p.num_compete_iter << std::endl
      << "Selection Rate, " << p.selection_rate << std::endl
      << "Mutation Rate, " << p.mutation_rate << std::endl
      << "Fitness Weights";
   for (auto w : p.weights)
      os << ", " << w;
   os << std::endl;
   return os;
}

std::istream& operator>>(std::istream & is, Simulation::Parameters& p)
{
   Simulation::Parameters tmp;
   is >> tmp.pop_size
      >> tmp.num_genetic_iter
      >> tmp.num_compete_iter
      >> tmp.selection_rate
      >> tmp.mutation_rate;
      
   for (size_t i = 0; i < tmp.weights.size(); ++i)
      is >> tmp.weights[i];
   
   if (is) {
      //ensure validity of ranges
      assert(tmp.pop_size > 2);
      assert(tmp.num_genetic_iter >= 1);
      assert(tmp.num_compete_iter >= 1);
      assert((tmp.selection_rate > 0.0) && (tmp.selection_rate <= 1.0));
      assert((tmp.mutation_rate >= 0.0) && (tmp.mutation_rate <= 1.0));
      
      p = std::move(tmp);
   }
   return is;
}

int main()
{
   Simulation::Parameters params;
   while (std::cin >> params)
   {
      Simulation sim(params);
      sim.genetic();
   }
}