from decimal import Decimal

Decimal()


class Scale(Decimal):

    def __new__(cls, factor, name):
        instance = super().__new__(cls, factor)
        instance.name = name
        return instance

    def __repr__(self):
        return f"Scale: {self.name} | factor: {self}"


SCALE_1_TO_1000 = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("1000"), "1:1000 Scale"
)
SCALE_1_TO_500 = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("500"), "1:500 Scale"
)
SCALE_1_TO_200 = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("200"), "1:200 Scale"
)
SCALE_1_TO_100 = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("100"), "1:100 Scale"
)
SCALE_1_TO_50 = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("50"), "1:50 Scale"
)
SCALE_1_TO_20 = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("20"), "1:20 Scale"
)
SCALE_1_TO_10 = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("10"), "1:10 Scale"
)

SCALE_1_TO_1000_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("1000") / Decimal("1000"),
    "1:1000 Scale, in meters",
)
SCALE_1_TO_500_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("500") / Decimal("1000"),
    "1:500 Scale, in meters",
)
SCALE_1_TO_200_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("200") / Decimal("1000"),
    "1:200 Scale, in meters",
)
SCALE_1_TO_100_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("100") / Decimal("1000"),
    "1:100 Scale, in meters",
)
SCALE_1_TO_50_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("50") / Decimal("1000"),
    "1:50 Scale, in meters",
)
SCALE_1_TO_20_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("20") / Decimal("1000"),
    "1:20 Scale, in meters",
)
SCALE_1_TO_10_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("25.4") * Decimal("10") / Decimal("1000"),
    "1:10 Scale, in meters",
)

SCALE_32ND_INCH = Scale(Decimal("1") / Decimal("72") * Decimal("32"), '1/32" = 1\'-0"')
SCALE_16TH_INCH = Scale(Decimal("1") / Decimal("72") * Decimal("16"), '1/16" = 1\'-0"')
SCALE_EIGHTH_INCH = Scale(Decimal("1") / Decimal("72") * Decimal("8"), '1/8" = 1\'-0"')
SCALE_QUARTER_INCH = Scale(Decimal("1") / Decimal("72") * Decimal("4"), '1/4" = 1\'-0"')
SCALE_HALF_INCH = Scale(Decimal("1") / Decimal("72") * Decimal("2"), '1/2" = 1\'-0"')
SCALE_ONE_INCH = Scale(Decimal("1") / Decimal("72"), '1" = 1\'-0"')

SCALE_3_32ND_INCH = Scale(
    Decimal("1") / Decimal("72") * Decimal("32") / Decimal("3"), '3/32" = 1\'-0"'
)
SCALE_3_16TH_INCH = Scale(
    Decimal("1") / Decimal("72") * Decimal("16") / Decimal("3"), '3/16" = 1\'-0"'
)
SCALE_3_EIGHTH_INCH = Scale(
    Decimal("1") / Decimal("72") * Decimal("8") / Decimal("3"), '3/8" = 1\'-0"'
)
SCALE_3_QUARTER_INCH = Scale(
    Decimal("1") / Decimal("72") * Decimal("4") / Decimal("3"), '3/4" = 1\'-0"'
)
SCALE_3_HALF_INCH = Scale(
    Decimal("1") / Decimal("72") * Decimal("2") / Decimal("3"), '3/2" = 1\'-0"'
)

SCALE_32ND_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("32") * Decimal("0.3048"),
    '1/32" = 1\'-0" IN METERS',
)
SCALE_16TH_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("16") * Decimal("0.3048"),
    '1/16" = 1\'-0" IN METERS',
)
SCALE_EIGHTH_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("8") * Decimal("0.3048"),
    '1/8" = 1\'-0" IN METERS',
)
SCALE_QUARTER_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("4") * Decimal("0.3048"),
    '1/4" = 1\'-0" IN METERS',
)
SCALE_HALF_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("2") * Decimal("0.3048"),
    '1/2" = 1\'-0" IN METERS',
)
SCALE_ONE_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("0.3048"), '1" = 1\'-0"'
)

SCALE_3_32ND_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("32") / Decimal("3") * Decimal("0.3048"),
    '3/32" = 1\'-0" IN METERS',
)
SCALE_3_16TH_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("16") / Decimal("3") * Decimal("0.3048"),
    '3/16" = 1\'-0" IN METERS',
)
SCALE_3_EIGHTH_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("8") / Decimal("3") * Decimal("0.3048"),
    '3/8" = 1\'-0" IN METERS',
)
SCALE_3_QUARTER_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("4") / Decimal("3") * Decimal("0.3048"),
    '3/4" = 1\'-0" IN METERS',
)
SCALE_3_HALF_INCH_M = Scale(
    Decimal("1") / Decimal("72") * Decimal("2") / Decimal("3") * Decimal("0.3048"),
    '3/2" = 1\'-0" IN METERS',
)
