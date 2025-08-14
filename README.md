# QRW_Markov_Walk

This repositorty contrains the behavoural data from 34 participants doing the human-AI teaming task outlined in:
https://www.mdpi.com/1099-4300/25/9/1362


### Data

On a trial-by-trail basis,
1. A crosshair appears for ~ 0.6 seconds - there is a bit of jutter (0.5 seconds approx).
2. The human particpant looks at a face, determines if it's real or fake, and then presses the corresponding button (left or right). The face is grayed out, and a blue or yellow avatar (of the human) is shown.
3. The 'AI' has in the meantime also made an assessment if it's a real or fake face, and presents its findings to the participant 300 msec after the human participant pressed. The AI avatar appears (yellow or blue).
4. For 1 second both icons remain on the screen before moving on to the next trial.

We have for each of these trials, the trial type as a code:
- 3 = "No Response"
- 4 = "HRAIR", i.e. trials in which the human pressed real and the AI also said the face it's real. (match trial)
- 5 = "HRAIF", i.e. human-real, AI-fake; a mismatch trial
- 6 = "HFAIF"
- 7 = "HFAIR"

For each trial, we also have the time at which the AI presented the feedback (so we are able to calculate the ITI from that).

For each trial, we finally provide the ratings. After every 28 trials, we give a VAS (visual analog scale) slider on the screen from 0 to 100%, with the question:
'How reliable do you find the AI'?

In total, the experiment has 28 Pratice trials (1 practice block), followed by 560 epxerimental trials (over 20 blocks - each block containing 28 trials).
There are 22 ratings in total (the rating before the practice block; rating after the practice block, and then 20 ratings after each experimental block.

### Scripts

This repository provides the Matlab scripts used to model the reliability ratings as a function of the response times and the event types, using the Quantum Random Walk approach and the Markov approach.
For the modeling paper, see:
https://arxiv.org/abs/2504.13918
This purely behavoural paper has been accepted in Philosophical Transaction A.


