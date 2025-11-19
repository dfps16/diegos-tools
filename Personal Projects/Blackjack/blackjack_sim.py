import random

class Card:
    suit = ['Hearts', 'Diamonds', 'Clubs', 'Spades']
    rank = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
            'Jack', 'Queen', 'King', 'Ace']

    def __init__(self, suit, rank):
        self.suit = suit
        self.rank = rank

        def __str__(self):
            return f"{self.suit} of {self.rank}"

class Deck:
    def __init__(self):
        self.cards=[Card(rank, suit), ]