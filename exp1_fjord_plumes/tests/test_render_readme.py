from exp1_fjord_plumes.scripts.render_readme import typographic_minus


def test_typographic_minus_signs_only():
    assert typographic_minus("-1.5") == "−1.5"
    assert typographic_minus("(-2.5") == "(−2.5"
    assert typographic_minus("0.50-0.61") == "0.50-0.61"
    assert typographic_minus("a_line2d-3") == "a_line2d-3"
    assert typographic_minus("pending") == "pending"
