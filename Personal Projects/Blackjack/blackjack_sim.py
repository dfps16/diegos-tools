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
    def __init__(self, num_decks=1):
        self.cards = [Card(rank, suit)
                      for i in range(num_decks)
                      for suit in Card.suit
                      for rank in Card.rank]
        random.shuffle(self.cards)

    def deal(self):
        return self.cards.pop()

    def __len__(self):
        return len(self.cards)


class Hand:
    def __init__(self):
        self.cards = []

    def add_card(self, card):
        self.cards.append(card)

    def total(self):
        value = 0
        aces = 0
        for card in self.cards:
            if card.rank in ['Jack', 'Queen', 'King']:
                value += 10
            elif card.rank == 'Ace':
                aces += 1
                value += 11
            else:
                value += int(card.rank)
            while value > 21 and aces > 0:
                value -= 10
                aces -= 1
        return value
    def __str__(self):
        return ', '.join(str(card) for card in self.cards)


class BlackjackGame:
    def __init__(self, num_decks=1):
        self.deck = Deck(num_decks)
        self.player_hand = Hand()
        self.dealer_hand = Hand()

    def start(self):
        print("Welcome to Blackjack!\n")
        playing = True

        # Initial deal
        for i in range(2):
            self.player_hand.add_card(self.deck.deal())
            self.dealer_hand.add_card(self.deck.deal())

        while playing:
            print(f"Your Hand: {self.player_hand} (Total: {self.player_hand.total()})")
            time.sleep(1)
            print(f"Dealer shows: {self.dealer_hand.cards[0]}\n")
            time.sleep(1)
            if self.player_hand.total() == 21:
                print("Blackjack! You win!")
                return
            elif self.player_hand.total() > 21:
                print("Bust! Dealer wins.")
                return

            time.sleep(1)
            action = input("Do you want to [h]it or [s]tand? ").lower()
            if action == 'h':
                self.player_hand.add_card(self.deck.deal())
            elif action == 's':
                playing = False
        time.sleep(1)
        print(f"Dealer's Hand: {self.dealer_hand} (Total: {self.dealer_hand.total()})")
        time.sleep(1)
        while self.dealer_hand.total() < 17:
            self.dealer_hand.add_card(self.deck.deal())
            print(f"Dealer hits: {self.dealer_hand} (Total: {self.dealer_hand.total()})")

        # Determine outcome
        player_total = self.player_hand.total()
        dealer_total = self.dealer_hand.total()
        if dealer_total > 21 or player_total > dealer_total:
            print("You win!")
        elif dealer_total == player_total:
            print("Push (tie)!")
        else:
            print("Dealer wins.")


if __name__ == "__main__":
    num_decks = int(input("How many decks? "))
    BlackjackGame(num_decks=num_decks).start()
