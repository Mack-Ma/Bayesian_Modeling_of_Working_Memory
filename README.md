# Bayesian Modeling of Working Memory (BMW) Toolbox   
### Ver. 2019

![](https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory/blob/master/BMW_icon.png)

### Prologue
Given the fast evolution of working memory computational models and the methodology of modeling these years, a new toolbox that could estimate, assess and comprehensively compare the new models using the state-of-the-art methods should be fruitful. Bayesian Modeling of Working Memory (BMW) toolbox is working on this issue by establishing a pipeline that integrates model definition, model estimation & model comparison. All users could either finish all the job by GUI or by recruiting sample scripts.

### Pipeline

- Fit Models by Maximum Likelihood (MLE)/ Maximum A Posteriori (MAP)
- Available Methods for Model Comparison
	- Log Likelihood (LLH)
	- Akaike Information Criterion (AIC)
	- Modified Akaike Information Criterion (AICc)
	- Bayesian Information Criterion (BIC)
	- Deviance Information Criterion (DIC/DIC*)  
	- Watanabe-Akaike Information Criterion (WAIC)  
	- Log Model Evidence/Marginal Likelihood (LME)
	- Second-Level Random-Effect Bayesian Model Selection (RFX-BMS)
- Available Models 
	- Continuous Recall Models
		- Item Limit
		- Standard Mixture
		- Slots-plus-Averaging
		- Equal Precision
		- Variable Precision
		- Variable Precision with Capacity
		- Categorical Slots-plus-Averaging   
		(Categorical encoding will appear as model variant in the future release)
		- Categorical Variable Precision
		- Categorical Variable Precision with Capacity
		- Model Variants: Bias/Swap/Bias Fluctuation/Precision Fluctuation
	- Change Detection Models
		- Fixed-Capacity (Single-Probe)  
		- Fixed-Capacity (Center-Probe)
		- Fixed-Capacity (Whole-Display)
		- Signal Detection
		- Model Variants: Lapse/Ensemble Encoding
	- Custom Models
		- Instance Model for Category Learning
		- Von Mises Distribution
- Available Optimization Algorithms  
	- Default: Differential Evolution Monte Carlo Markov Chain (DE-MCMC)  
  See "Requirements" for other options
  
### Requirements

- Matlab (Best if >= 9.0/2016a, not sure if the lower version works)
- Optimization Algorithm
  - __Default:__ _**Differential Evolution Monte Carlo Markov Chain (DE-MCMC)**_, built-in function in BMW toolbox
  - _**Adaptive Metropolis-Hastings Monte Carlo Markov Chain (MH-MCMC)**_, built-in function in BMW toolbox
  - _**fmincon**_ (sqp/interior point/active set), Optimization Toolbox (built-in toolbox in Matlab >= 9.0)
  - _**Genetic Algorithm**_, Global Optimization Toolbox (built-in toolbox in Matlab >= 7.0)
  - _**Simulated Annealing**_, Global Optimization Toolbox (built-in toolbox in Matlab >= 7.0)
  - _**Mesh Adaptive Direct Search (MADS)**_, Global Optimization Toolbox (built-in toolbox in Matlab >= 7.0)
  - _**Bayesian Adaptive Direct Search (BADS)**_, BADS Toolbox (http://github.com/lacerbi/bads)

### Contributors  
**Programmed** by Ma, Tianye  
**Mentored** by Yixuan Dr. Ku (Sun Yat-Sen University)  
**Sincere thanks** to  
Sizhu Han,  
Ruyuan Dr. Zhang,  
Kuangshi Zhao,  
__for their invaluable help!__
  
### Epilogue
Computational modeling could act as a powerful tool especially in understanding abstract issues such as "Do we have discrete capacity in working memory storage?" and "How is our working memory system organized?". And we suppose that models of best quality should on the one hand take a good balance of accuracy and complexity while on the other hand offer useful outputs. That is to say, we take the well-renowned George & Draper's (1987) view,

> All models are wrong, but some are useful.

In order to make this toolbox a convenient tool and a good communicator of computational methods and working memory models for cognitive scientists, we'll always need feedback. So please feel free to start a Github issue or contact Ma, Tianye (mack_ma2018@outlook.com) for bug reports and any other kinds of advice.
