\ Model RTEcnomic_dispatch
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  14 Power_generation0 + 15 Power_generation1 + 30 Power_generation2
   + 40 Power_generation3 + 10 Power_generation4
Subject To
 TC_p4: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   - 0.3489894581346941 Power_generation2
   + 0.1595380380044316 Power_generation4 <= 152.534747868126
 TC_n4: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   - 0.3489894581346941 Power_generation2
   + 0.1595380380044316 Power_generation4 >= -647.465252131874
 TC_p4: 0.4375881286510441 Power_generation0
   + 0.4375881286510441 Power_generation1
   + 0.1894514201302626 Power_generation2
   + 0.3600101778652635 Power_generation4 <= 1.0000001343382797e+09
 TC_n4: 0.4375881286510441 Power_generation0
   + 0.4375881286510441 Power_generation1
   + 0.1894514201302626 Power_generation2
   + 0.3600101778652635 Power_generation4 >= -9.9999986566172028e+08
 TC_p4: 0.3684952662324581 Power_generation0
   + 0.3684952662324581 Power_generation1
   + 0.1595380380044316 Power_generation2
   - 0.5195482158696954 Power_generation4 <= 1.0000001131269724e+09
 TC_n4: 0.3684952662324581 Power_generation0
   + 0.3684952662324581 Power_generation1
   + 0.1595380380044316 Power_generation2
   - 0.5195482158696954 Power_generation4 >= -9.9999988687302756e+08
 TC_p4: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   - 0.3489894581346944 Power_generation2
   + 0.1595380380044316 Power_generation4 <= 1.0000000525347478e+09
 TC_n4: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   - 0.3489894581346944 Power_generation2
   + 0.1595380380044316 Power_generation4 >= -9.9999994746525216e+08
 TC_p4: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   + 0.6510105418653058 Power_generation2
   + 0.1595380380044317 Power_generation4 <= 1.0000003525347478e+09
 TC_n4: 0.1939166051164977 Power_generation0
   + 0.1939166051164977 Power_generation1
   + 0.6510105418653058 Power_generation2
   + 0.1595380380044317 Power_generation4 >= -9.9999964746525216e+08
 TC_p4: - 0.3684952662324582 Power_generation0
   - 0.3684952662324582 Power_generation1
   - 0.1595380380044316 Power_generation2
   - 0.480451784130305 Power_generation4 <= 126.8730275968576
 TC_n4: - 0.3684952662324582 Power_generation0
   - 0.3684952662324582 Power_generation1
   - 0.1595380380044316 Power_generation2
   - 0.480451784130305 Power_generation4 >= -353.1269724031424
 R12: Power_generation0 + Power_generation1 + Power_generation2
   + Power_generation3 + Power_generation4 = 1000
Bounds
 Power_generation0 <= 40
 Power_generation1 <= 170
 Power_generation2 <= 520
 Power_generation3 <= 200
 Power_generation4 <= 600
End
