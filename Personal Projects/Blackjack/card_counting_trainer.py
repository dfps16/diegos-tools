import random
import time

class Card:
    suit = ['Hearts', 'Diamonds', 'Clubs', 'Spades']
    rank = ['2', '3', '4', '5', '6', '7', '8', '9', '10',
            'Jack', 'Queen', 'King', 'Ace']

    def __init__(self, rank, suit):
        self.suit = suit
        self.rank = rank

    def __str__(self):
        return f"{self.rank} of {self.suit}"


class Deck:
    def __init__(self):
        self.cards = [Card(rank, suit) for suit in Card.suit for rank in
                      Card.rank]
        random.shuffle(self.cards)

    def deal(self):
        return self.cards.pop()


class Trainer:
    def __init__(self):
        self.deck = Deck()

    def train(self):
        print("Beginning Training\n")
        playing = True
        while playing:
            card_no = input("Number of cards to deal: ")
            wait = input("Time (s) between cards: ")
            rolling_count = 0
            for i in range(int(card_no)):
                card = self.deck.deal()
                if card.rank in ['2', '3', '4', '5', '6']:
                    rolling_count += 1
                elif card.rank in ['10', 'Jack', 'Queen', 'King', 'Ace']:
                    rolling_count += -1
                print(card)
                time.sleep(float(wait))

            guess = int(input("What's the count?: "))
            if rolling_count == guess:
                print(f'You are correct! The count is {rolling_count}.')
            else:
                print(f'You are incorrect! The count is {rolling_count}.')
            game_end = True
            while game_end:
                encore = input("Do you want to train again? (y/n): ")
                if encore == 'n':
                    playing = False
                    game_end = False
                elif encore == 'y':
                    playing = True
                    game_end = False
                else:
                    print("That's not a valid option.")



if __name__ == "__main__":
    Trainer().train()