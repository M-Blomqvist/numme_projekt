function [U_star,max_star]= find_period(wanted_period,guess1,guess2, start, h,stop, interpol_steps_per_step)
    diff = 100;
    while diff > 1e-10 
        fprintf('\n guesses = %d & %d, diff = %d\n', guess1,guess2,diff);
        [~,~,~,y_max,i_period] = interpol(guess1, start, h, stop,interpol_steps_per_step);
        fprintf('\n interpol %d\n', guess2);
        [~,~,~,~,i_period_h] = interpol(guess2, start, h, stop, interpol_steps_per_step);
        fprintf('\n i_period = %d, i_period_h = %d\n', i_period,i_period_h);
        guess = guess1-(i_period*(guess1-guess2))/(i_period-i_period_h);
        guess1 = guess2;
        guess2 = guess;
        diff = abs(i_period - wanted_period);
    end
    U_star = guess;
    max_star = y_max;
end