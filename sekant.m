function root = sekant(f,guess1,guess2,certainty)
    diff = 100;
    fprintf('\n Sekant: guesses = %d & %d \n', guess1,guess2);
    while diff > certainty
        y1 = f(guess1);
        y2 = f(guess2);
        guess = guess1-(y1*(guess1-guess2))/(y1-y2);
        guess1 = guess2;
        guess2 = guess;
        diff = abs(guess1-guess2);
        %fprintf('\n guesses = %d & %d, diff: %d \n', guess1,guess2, diff);
    end
    root = guess1;
end