function root = sekant(f,guess1,guess2,certainty)
    diff = 100;
    
    guess_1 = guess1;
    guess_2 = guess2;
    %bryt då skillanden i guesses är under certainty
    while diff > certainty
        y1 = f(guess_1);
        y2 = f(guess_2);
        %sekant metoden
        guess = guess_1-(y1*(guess_1-guess_2))/(y1-y2);
        
        guess_1 = guess_2;
        guess_2 = guess;
        diff = abs(guess_1-guess_2);
    end
    root = guess_1;
end