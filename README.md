# HMMETS
HMM for Evolution of Targeting Signal (HMMETS) is an estimation algorithm for targeting signal of proteins.
Basic idea came from phylogenetic profiling for gene gain/loss, and we expanded this to targeting signal.
One motivation of this algorithm is finding widely conserved targeting signal in a gene cluster.
Concept of Homology is hard to be applied in this task, since targeting signal is usually not conserved in a promary structure level.
However, targeting signal of orthologs is basically conserved, inferring conservation not in primary sequence but in feature space.
Assuming there is no re-entrant into feature space of targeting signal from outside, we can simply model evolution of targeting signal with some prediction algotirhms.
At present, we applied this for estimation of presequence evolution with MitoFates program.
## Installation
This projects is coded in ruby, so you don't need to compile.
Download and enjoy!
Check Gemfile for dependency.

Dependency:
1. bioruby
2. prawn
3. GSL (required ruby binding)
4. RubyInline
5. parallel

## Usage
Writing
## Contributing
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D
## Credits
Yoshinori Fukasawa
