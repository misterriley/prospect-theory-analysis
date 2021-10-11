functions
{
  real get_weight(real prob, real delta)
  {
    real numerator;
    real denominator;
    
    if(prob == 0)
    {
      return 0;
    }
      
    if(prob == 1)
    {
      return 1;
    }
    
    if(delta == 1)
    {
      return prob;
    }
    
    if(delta == 0)
    {
       return .5;
    }
    
    numerator = prob^delta;
    denominator = (numerator + (1-prob)^delta)^(1/delta);
    
    return numerator/denominator;
  }
  
  real get_utility(real gain, real alpha, real lambda)
  {
    if(gain == 0)
    {
      return 0;
    }
    
    if(gain < 0)
    {
      return -1 * lambda * (-1 * gain)^alpha;
    }
    
    return gain^alpha;
  }
  
  real get_pt_value(real[,] codes, int code_index, real alpha, real delta, real lambda)
  {
    real gain1; 
    real prob1; 
    real gain2; 
    real prob2; 
    real value1; 
    real value2;
    
    gain1 = codes[code_index,3];
    prob1 = codes[code_index,4];
    gain2 = codes[code_index,5];
    prob2 = codes[code_index,6];
    
    value1 = get_utility(gain1, alpha, lambda) * get_weight(prob1, delta);
    value2 = get_utility(gain2, alpha, lambda) * get_weight(prob2, delta);
    
    return value1 + value2;
  }
}
data
{
  int<lower=1> N; //Number of subjects
  int<lower=0,upper=2> choices[N,16];
  real codes[32,6];
}
parameters
{
  real <lower=0> alpha_sd;
  real <lower=0, upper=1> alphas[N];
  
  //real <lower=0.01, upper=pi()/2> c_mean;
  real <lower=0> c_sd;
  real <lower=0.0001, upper=pi()/2> cs[N];

  real <lower=0> delta_sd;
  real <lower=0,upper=1> deltas[N];
  
  real <lower=0>lambda_sd;
  real <lower=0.0001>lambdas[N];
}
model
{
  real inv_temp;
  real pt_value1;
  real pt_value2;
  
  int choice;
  int option1Index;
  int option2Index;
  int outcome;
  
  //hyperparameters have densities taken from hyperpriors
  alpha_sd ~ cauchy(0, 1);
  c_sd ~ cauchy(0, 1);
  delta_sd ~ cauchy(0, 1);
  lambda_sd ~ cauchy(0, 1);
  
  //iterate over subjects
  for(s in 1:N)
  {
    //parameters have normal densities from priors, median estimates taken from Kahneman & Tversky (1992)
    alphas[s] ~ normal(.88, alpha_sd);
    cs[s] ~ normal(pi()/4, c_sd);
    deltas[s] ~ normal(.65, delta_sd);
    lambdas[s] ~ normal(2.25, lambda_sd);
    
    inv_temp = atan(cs[s]);
    
    //iterate over frames
    for(f in 1:16)
    {
      choice = choices[s,f];
      
      //choice == 0 is a recoding of NA in the original file - skip
      if(choice != 0)
      {
        option1Index = 2*f - 1;
        option2Index = 2*f;
        
        pt_value1 = get_pt_value(codes, option1Index, alphas[s], deltas[s], lambdas[s]);
        pt_value2 = get_pt_value(codes, option2Index, alphas[s], deltas[s], lambdas[s]);
        
        //Choice is "1" or "2", so we subtract 1 to make it 0 or 1.   
        outcome = choice - 1;
        
        //If pt_value2 > pt_value1, then the participant should have chosen 2.  Thus, outcome should be 1,
        //and 1 ~ bernoulli_logit(inv_temp * (pt_value2 - pt_value1)) should have a high posterior.
        //Conversely, if pt_value2 < pt_value1, then 0 ~ bernoulli_logit(...) should have a high posterior.  
        outcome ~ bernoulli_logit(inv_temp * (pt_value2 - pt_value1));
      }
    }
  }
}
