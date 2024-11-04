# flyout_prob
Baseball flyout probability model + trajectory calculator for larger optimal fielder positioning paper

- My contribution as second author to "Optimizing Baseball Fielder Positioning with Consideration for Adaptable Hitters" paper found here, https://www.sloansportsconference.com/research-papers/optimizing-baseball-fielder-positioning-with-consideration-for-adaptable-hitters
- Calculated flight path of batted balls in R using batted ball features from baseballR, basing physics calculations off Alan M. Nathan's trajectory calculator
- Using calculated hang time and xy coordinates, created likelihood model using avergae fielder positioning, speed/acceleration, jump to get base flyout probabilities to create csv
- With csv, ran logistic regression models to find best combination of features
