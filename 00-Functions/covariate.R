covariate <- function(data, s, r, stp){
  
  ##############################################################################
  ##############################################################################
  if( stp < dat.gam$stp[1]){
    return(rep(0,8))
  } else {
    ind = max(which(data$stp <= stp))
  }
  
  # check if the transaction with same sender and receiver has already occurred
  prev <- which(data[1:(ind-1),"s1"]==s & data[1:(ind-1),"r1"]==r)
  
  # IF THE BANK S ALREADY BORROWED MONEY TO BANK R
  # then dynamics are analysed after that event
  # STARTING POINT SET AFTER THE LAST OF THESE EVENTS
  # (important choice from a dynamical pov)
  # inertia indicator is set equal to 1
  # time from the last inertial event is set equal to
  # the time of the latest transaction of this type
  if (length(prev)>0){
    inert.id <- 1
    exp.inert <- exp(-(stp-max(data[prev,"stp"])))
    prev <- prev[length(prev)] + 1
    
    # IF NOT
    # STARTING POINT SET EQUAL TO THE BEGINNING OF THE EVENT SEQUENCE
    # indicator and the exp.time covariates are set equal to 0
  } else {
    inert.id <- 0
    exp.inert <- 0
    prev <- 1
  }
  
  # CHECK IF THE BANK THAT IS RECEVING MONEY ALREADY RECEIVED MONEY
  # FROM SOME OTHER BANKS
  # Namely, check the transactions between starting point and the current time
  # with as receiver the receiver of current transaction
  t <- which(data[prev:ind,"r1"]==r)
  t <- prev + t - 1
  
  # STORE THOSE WHO SENT MONEY TO THE BANK THAT IS CURRENTLY BORROWING
  # the current sender is not included because the starting point
  # is set after the last event of this type
  l <- data[t,2]
  common.nodes <- unique(l)
  
  # INITIALIZE THE OTHER COVARIATES
  exp.triadic <- triadic.id <- 0
  exp.cyclic <- cyclic.id <- 0
  exp.triad <- triad.id <- 0
  exp.rec <- rec.id <- 0
  
  # IF AT LEAST ONE BANK (different from s) SENT MONEY TO THE BORROWING BANK
  if(length(common.nodes)!=0){
    
    # we inspect each of them
    for(cn in common.nodes){
      
      # if more than one event occurred we consider the last one
      stop_t <- data$stp[t[which(data[t,"s1"]==cn & data[t,"r1"] == r)]]
      stop_t <- max(stop_t)
      
      # we need to understand if this bank received money from the sender
      # in a time window
      # after the starting point
      # but before that the third part sent money to the receiver
      t2 <- which(data[prev:max(which(data$stp==stop_t)),"s1"]==s &
                    data[prev:max(which(data$stp==stop_t)),"r1"]==cn)
      
      # In case at least one event satisfies these conditions
      # the the indicator for triadic closure is set equal to 1
      
      # TO RECAP: 
      # triadic closure is set to 1 if:
      # the third part bank has already borrowed money to the receiver
      # the third part bank already received money by the sender
      
      if (length(t2)>0){
        triadic.id <- max(triadic.id, 1)
        exp.triadic <- max(exp.triadic, exp(-(stp-stop_t)))
      } else {
        t3 <- which(data[1:prev,"s1"]==s & data[1:prev,"r1"]==cn)
        # the event of the third part bank receiving money by the sender
        # is allowed to be recognized before the starting point 
        # but from the beginning of the analysis as well
        # the triadic closure time is computed by referring to the
        # time of the third part sending money to the current borrower
        if (length(t3)>0){
          triadic.id <- max(triadic.id, 1)
          exp.triadic <- max(exp.triadic, exp(-(stp-stop_t)))
        }
      }
    }
  }
  
  ##############################################################################
  ############################################################################## 
  
  # Once the starting point of the analysis has been defined, 
  
  # CHECK IF THE BANK THAT IS SENDING MONEY RECEIVED MONEY PREVIOUSLY
  # Namely, check the transactions between starting point and the current time
  # with as receiver the sender of current transaction
  t <- which(data[prev:ind,"r1"]==s)
  t <- prev + t - 1
  
  # STORE THOSE WHO SENT MONEY TO THE BANK THAT IS CURRENTLY BORROWING
  l <- data[t,2]
  common.nodes <- unique(l)
  
  # IF AT LEAST ONE BANK DID SENT MONEY
  if(length(common.nodes)!=0){
    
    # we inspect each of them
    for(cn in common.nodes){
      
      ##########################################################################
      # if this bank coincides with the receiver of the current event
      # then RECIPROCITY is there
      # the receiver has already borrowed money to the sender
      
      if(cn==r){
        
        # check if more than one reciprocal event has already occurred
        # if more than one store the latest time
        stop_t <- data$stp[t[which(data[t,"s1"]==cn & data[t,"r1"] == s)]]
        stop_t <- max(stop_t)
        rec.id <- max(rec.id,1)
        exp.rec <- exp(-(stp-stop_t))
        
        # if this bank does not coincide with the receiver of the event
        # then we need to check if there might be a cyclic closure
        
      } else {
        
        # in this case there is a third part bank that is involved
        # and that borrowed money to the bank that is now borrowing
        # we store the time of the last transaction 
        # in case more than one are there
        stop_t <- data$stp[t[which(data[t,"s1"]==cn & data[t,"r1"] == s)]]
        stop_t <- max(stop_t)
        
        # we need to understand if this bank received money
        # from the bank that is receiving money
        # in a time window
        # after the starting point
        # but before that the third part sent money to the sender
        t2 <- which(data[prev:max(which(data$stp<=stop_t)),"s1"]==r &
                      data[prev:max(which(data$stp<=stop_t)),"r1"]==cn)
        
        # In case at least one event satisfies these conditions
        # the the indicator for cyclic closure is set equal to 1
        
        # TO RECAP: 
        # cyclic closure is set to 1 if:
        # the third part bank has already borrowed money to the sender
        # the third part bank already received money by the receiver
        
        # The event might close more than one triangle 
        # The variable n.triad has the role
        # of considering the number of cycles that are actually closed
        if (length(t2)>0){
          cyclic.id <- max(cyclic.id, 1)
          exp.cyclic <- max(exp.cyclic, exp(-(stp-stop_t)))
        } else {
          t3 <- which(data[1:prev,"s1"]==r & data[1:prev,"r1"]==cn)
          # the event of the third part bank sending money to the receiver
          # is allowed to be recognized before the starting point 
          # but from the beginning of the analysis as well
          # the cyclic closure time is computed by referring to the
          # time of the third part sending money to the current lender
          if (length(t3)>0){
            cyclic.id <- max(cyclic.id, 1)
            exp.cyclic <- max(exp.cyclic, exp(-(stp-stop_t)))
          }
        }
      }
    }
    
    ############################################################################
    
    # if no bank already borrowed money to the sender bank
    # then the event is not closing either reciprocity or cyclic closure
  } else {
    rec.id <- 0
    exp.rec <- 0
    cyclic.id <- 0
    exp.cyclic <- 0
  }
  
  # -    Indicator for reciprocity;
  # -    Exp. of minus difference in time from last occurred reciprocal event;
  # -    Indicator for cyclic closure;
  # -    Exp. of minus difference in time from last occurred event of interest for cc;
  # -    Indicator for inertia;
  # -    Exp. of minus difference in time from last occurred inertial event;
  # -    Indicator for triadic closure;
  # -    Exp. of minus difference in time from last occurred event of interest for tc;
  
  return(c(rec.id, exp.rec, 
           cyclic.id, exp.cyclic, 
           inert.id, exp.inert, 
           triadic.id, exp.triadic))
}
