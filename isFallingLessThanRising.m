function value = isFallingLessThanRising(falling,rising)
if(length(falling) == 0)
  value = false;
elseif(length(rising) ==0)
    value = true;
else
  if(falling(1) < rising(1))
    value = true;
  else
    value = false;
  end  
end
end