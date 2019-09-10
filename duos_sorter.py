#!/usr/bin/env python3
import argparse

def find_duos(players_list):
    
    # Make a dictionay of each player which contains a dictionary of played with 
    # and vsed each of which contains a list of players they've played with or 
    # vsed respectively
    players_dict = {}
    for play in players_list:
        players_dict[play] = {"played_with" : [], "played_against" : []}

    #TODO add thing that removes spare players
    playing_this_round = players_list


    new_matchs = True
    round_number = 1
    while new_matchs:
        # Work out all possible pairs
        possible_pairs = []
        for p1i, player_1 in enumerate(playing_this_round[:-1]):
            for player_2 in playing_this_round[(p1i + 1):]:
                # If you imagine a matrix this is just the top right not including the diagonal
                possible_pairs.append([player_1, player_2])

        # Work out all possible wais the pairs could vs
        possible_matches = []
        for pair_1 in possible_pairs:
            for pair_2 in possible_pairs:
                #check for inpossible pairings
                if not (pair_1[0] in pair_2 or pair_1[1] in pair_2):
                    possible_matches.append([pair_1, pair_2])

        # Score each possible match by counting how many times the players have previous played against or played with their matches
        match_score_with = []
        match_score_against = []
        for match in possible_matches:
            temp_score_with = 0
            temp_score_against = 0
            pair_1, pair_2 = match
            py1, py2 = pair_1
            py3, py4 = pair_2
            # Check if played with
            if py1 in players_dict[py2]["played_with"]:
                temp_score_with += 1
            if py3 in players_dict[py4]["played_with"]:
                temp_score_with += 1

            # Check if played against
            if py1 in players_dict[py3]["played_against"]:
                temp_score_against += 1
            if py1 in players_dict[py4]["played_against"]:
                temp_score_against += 1
            if py2 in players_dict[py3]["played_against"]:
                temp_score_against += 1
            if py2 in players_dict[py4]["played_against"]:
                temp_score_against += 1
            #print("          Pair 1: {0}   {1}".format(py1, py2))
            #print("          Pair 2: {0}   {1}".format(py3, py4))
            #print("          score with {0}".format(temp_score_with))

            
            match_score_with.append(temp_score_with)
            match_score_against.append(temp_score_against)

        # Print the first match with the lowest score
        min_score_with = min(match_score_with)
        temp = []
        for mi, msa in enumerate(match_score_against):
            if match_score_with[mi] == min_score_with:
                temp.append(msa)
        min_score_against = min(temp)
        if min_score_with > 0:
            print("Out of possible matches")
            break
        for mi, match in enumerate(possible_matches):
            if match_score_with[mi] == min_score_with and \
               match_score_against[mi] == min_score_against:
                # Print the next match
                print("Round #{0}:".format(round_number))
                #for ci, court in enumerate(match):
                pair_1, pair_2 = match
                py1, py2 = pair_1
                py3, py4 = pair_2
                print("          Pair 1: {0}   {1}".format(py1, py2))
                print("          Pair 2: {0}   {1}".format(py3, py4))

                # Update the dicts
                round_number += 1
                players_dict[py1]["played_with"].append(py2)
                players_dict[py2]["played_with"].append(py1)
                players_dict[py3]["played_with"].append(py4)
                players_dict[py4]["played_with"].append(py3)

                players_dict[py1]["played_against"].append(py3)
                players_dict[py1]["played_against"].append(py4)
                players_dict[py2]["played_against"].append(py3)
                players_dict[py2]["played_against"].append(py4)
                players_dict[py3]["played_against"].append(py1)
                players_dict[py3]["played_against"].append(py2)
                players_dict[py4]["played_against"].append(py1)
                players_dict[py4]["played_against"].append(py2)

                #print(players_dict)
                #if round_number == 3:
                #    exit()
                break

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
            Sorts a bunch of players into duos with the least number of people vsed and played with.
            """)
    parser.add_argument("-p", "--players", nargs='*', type=str,
                        help="Space seperated list of names. Eg: Nick Hannah Keegan")
    args=parser.parse_args()
    
    find_duos(args.players)
