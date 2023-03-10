
test_that('helper_wZ',
          {
            roots <- wZ(p = p, weight=weight, is.onetail = is.onetail, eff.sign = eff.sign)
            expect_lte(roots$p,1)
          })
test_that('helper_wFisher',
          {
            roots <- wZ(p = p, weight=weight, is.onetail = is.onetail, eff.sign = eff.sign)
            expect_lte(roots$p,1)
          })
test_that('helper_lancaster',
          {
            roots <- wZ(p = p, weight=weight, is.onetail = is.onetail, eff.sign = eff.sign)
            expect_lte(roots$p,1)
          })
