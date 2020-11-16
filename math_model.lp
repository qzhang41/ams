\ Model RTEcnomic_dispatch
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  14 Power_generation0 + 15 Power_generation1 + 30 Power_generation2
   + 40 Power_generation3 + 10 Power_generation4
Subject To
 TC_p0: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   - 0.3489894581346941 Power_generation2
   + 0.1595380380044316 Power_generation4 <= 294.8272678439536
 TC_n0: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   - 0.3489894581346941 Power_generation2
   + 0.1595380380044316 Power_generation4 >= -505.1727321560464
 TC_p1: 0.4375881286510441 Power_generation0
   + 0.4375881286510441 Power_generation1
   + 0.1894514201302626 Power_generation2
   + 0.3600101778652635 Power_generation4 <= 1.0000000570937688e+09
 TC_n1: 0.4375881286510441 Power_generation0
   + 0.4375881286510441 Power_generation1
   + 0.1894514201302626 Power_generation2
   + 0.3600101778652635 Power_generation4 >= -9.9999994290623116e+08
 TC_p2: 0.3684952662324581 Power_generation0
   + 0.3684952662324581 Power_generation1
   + 0.1595380380044316 Power_generation2
   - 0.5195482158696954 Power_generation4 <= 1.0000000480789633e+09
 TC_n2: 0.3684952662324581 Power_generation0
   + 0.3684952662324581 Power_generation1
   + 0.1595380380044316 Power_generation2
   - 0.5195482158696954 Power_generation4 >= -9.9999995192103672e+08
 TC_p3: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   - 0.3489894581346944 Power_generation2
   + 0.1595380380044316 Power_generation4 <= 1.0000000223272679e+09
 TC_n3: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   - 0.3489894581346944 Power_generation2
   + 0.1595380380044316 Power_generation4 >= -9.9999997767273211e+08
 TC_p4: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   + 0.6510105418653058 Power_generation2
   + 0.1595380380044317 Power_generation4 <= 1.0000001498272679e+09
 TC_n4: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   + 0.6510105418653058 Power_generation2
   + 0.1595380380044317 Power_generation4 >= -9.9999985017273211e+08
 TC_p5: - 0.3684952662324582 Power_generation0
   - 0.3684952662324582 Power_generation1
   - 0.1595380380044316 Power_generation2
   - 0.480451784130305 Power_generation4 <= 191.9210367286645
 TC_n5: - 0.3684952662324582 Power_generation0
   - 0.3684952662324582 Power_generation1
   - 0.1595380380044316 Power_generation2
   - 0.480451784130305 Power_generation4 >= -288.0789632713355
 balance: Power_generation0 + Power_generation1 + Power_generation2
   + Power_generation3 + Power_generation4 = 425
Bounds
 Power_generation0 <= 40
 Power_generation1 <= 170
 Power_generation2 <= 520
 Power_generation3 <= 200
 Power_generation4 <= 600
End
