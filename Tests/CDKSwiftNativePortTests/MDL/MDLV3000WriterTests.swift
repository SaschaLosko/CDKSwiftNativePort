import Foundation
#if canImport(CoreGraphics)
import CoreGraphics
#endif
import XCTest
@testable import CDKSwiftNativePort

final class MDLV3000WriterTests: XCTestCase {

    func testWritesBasicV3000Molfile() throws {
        let molecule = Molecule(
            name: "BasicV3000",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.2, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("999 V3000"))
        XCTAssertTrue(text.contains("M  V30 BEGIN ATOM"))
        XCTAssertTrue(text.contains("M  V30 1 C 0 0 0 0"))
        XCTAssertTrue(text.contains("M  V30 2 O 1.2 0 0 0"))
        XCTAssertTrue(text.contains("M  V30 1 1 1 2"))
    }

    func testOutputsValencyWhenNeeded() throws {
        let molecule = Molecule(
            name: "Valency",
            atoms: [
                Atom(id: 1, element: "Na", position: .zero, explicitHydrogenCount: 0),
                Atom(id: 2, element: "Na", position: CGPoint(x: 1.2, y: 0), explicitHydrogenCount: 1)
            ],
            bonds: []
        )

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  V30 1 Na 0 0 0 0 VAL=-1"))
        XCTAssertTrue(text.contains("M  V30 2 Na 1.2 0 0 0"))
        XCTAssertFalse(text.contains("M  V30 2 Na 1.2 0 0 0 VAL="))
    }

    func testWritesMassRadicalAndPushesHydrogensToBack() throws {
        let molecule = Molecule(
            name: "MassAndHydrogenOrder",
            atoms: [
                Atom(id: 1, element: "H", position: CGPoint(x: 0.0, y: 0.0), isotopeMassNumber: 2),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.2, y: 0.0), explicitHydrogenCount: 3, radical: 2)
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single)
            ]
        )

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  V30 1 C 1.2 0 0 0 RAD=2"))
        XCTAssertTrue(text.contains("M  V30 2 H 0 0 0 0 MASS=2"))
    }

    func testWritesFormalChargeAndReversedWedgeEndpoints() throws {
        let molecule = Molecule(
            name: "ChargeAndWedge",
            atoms: [
                Atom(id: 1, element: "O", position: .zero, charge: -1),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.2, y: 0.0), explicitHydrogenCount: 3),
                Atom(id: 3, element: "N", position: CGPoint(x: 2.4, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, stereo: .upReversed),
                Bond(id: 2, a1: 2, a2: 3, order: .single, stereo: .downReversed)
            ]
        )

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  V30 1 O 0 0 0 0 CHG=-1"))
        XCTAssertTrue(text.contains("M  V30 1 1 2 1 CFG=1"))
        XCTAssertTrue(text.contains("M  V30 2 1 3 2 CFG=3"))
    }

    func testWritesAtomAttributesAndQueryTypes() throws {
        let molecule = Molecule(
            name: "Attributes",
            atoms: [
                Atom(id: 1,
                     element: "C",
                     position: CGPoint(x: 0.0, y: 0.0),
                     charge: -1,
                     isotopeMassNumber: 13,
                     radical: 2),
                Atom(id: 2,
                     element: "L",
                     position: CGPoint(x: 1.4, y: 0.0),
                     explicitHydrogenCount: 0,
                     atomList: ["N", "O"],
                     atomListIsNegated: true),
                Atom(id: 3,
                     element: "R1",
                     position: CGPoint(x: 2.8, y: 0.0),
                     rGroupLabel: 1,
                     attachmentPoint: 1)
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .single)
            ]
        )

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("CHG=-1"))
        XCTAssertTrue(text.contains("MASS=13"))
        XCTAssertTrue(text.contains("RAD=2"))
        XCTAssertTrue(text.contains("ATOMLIST=(2 N O)"))
        XCTAssertTrue(text.contains("NOT=TRUE"))
        XCTAssertTrue(text.contains("HCOUNT=0"))
        XCTAssertTrue(text.contains("RGROUPS=(1 1)"))
        XCTAssertTrue(text.contains("ATTCHPT=1"))
    }

    func testWritesQueryBondTypesAndBondTopology() throws {
        let molecule = Molecule(
            name: "BondQueries",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.2, y: 0.0)),
                Atom(id: 3, element: "C", position: CGPoint(x: 2.4, y: 0.0)),
                Atom(id: 4, element: "C", position: CGPoint(x: 3.6, y: 0.0)),
                Atom(id: 5, element: "C", position: CGPoint(x: 4.8, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, stereo: .up),
                Bond(id: 2, a1: 2, a2: 3, order: .single, stereo: .either),
                Bond(id: 3, a1: 3, a2: 4, order: .single, queryType: .singleOrDouble, topology: .ring),
                Bond(id: 4, a1: 4, a2: 5, order: .single, queryType: .any, topology: .chain)
            ]
        )

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  V30 1 1 1 2 CFG=1"))
        XCTAssertTrue(text.contains("M  V30 2 1 2 3 CFG=2"))
        XCTAssertTrue(text.contains("M  V30 3 5 3 4 TOPO=1"))
        XCTAssertTrue(text.contains("M  V30 4 8 4 5 TOPO=2"))
    }

    func testWritesAromaticBondTypeFour() throws {
        let atoms = (1...6).map { index in
            Atom(id: index, element: "C", position: CGPoint(x: Double(index), y: 0.0))
        }
        let bonds = [
            Bond(id: 1, a1: 1, a2: 2, order: .aromatic),
            Bond(id: 2, a1: 2, a2: 3, order: .aromatic),
            Bond(id: 3, a1: 3, a2: 4, order: .aromatic),
            Bond(id: 4, a1: 4, a2: 5, order: .aromatic),
            Bond(id: 5, a1: 5, a2: 6, order: .aromatic),
            Bond(id: 6, a1: 6, a2: 1, order: .aromatic)
        ]
        let text = try CDKMDLV3000Writer.write(Molecule(name: "Aromatic", atoms: atoms, bonds: bonds))

        XCTAssertTrue(text.contains("M  V30 1 4 1 2"))
        XCTAssertTrue(text.contains("M  V30 6 4 6 1"))
    }

    func testSelectionWrittenAsHILITECollection() throws {
        var molecule = Molecule(
            name: "Highlighted",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 1.5000, y: 0.0000)),
                Atom(id: 2, element: "C", position: CGPoint(x: 0.7500, y: 1.2990)),
                Atom(id: 3, element: "C", position: CGPoint(x: -0.7500, y: 1.2990)),
                Atom(id: 4, element: "C", position: CGPoint(x: -1.5000, y: 0.0000)),
                Atom(id: 5, element: "C", position: CGPoint(x: -0.7500, y: -1.2990)),
                Atom(id: 6, element: "C", position: CGPoint(x: 0.7500, y: -1.2990))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .double),
                Bond(id: 3, a1: 3, a2: 4, order: .single),
                Bond(id: 4, a1: 4, a2: 5, order: .double),
                Bond(id: 5, a1: 5, a2: 6, order: .single),
                Bond(id: 6, a1: 6, a2: 1, order: .double)
            ]
        )
        molecule.highlightedAtomIDs = [2, 3, 5, 4]
        molecule.highlightedBondIDs = [2]

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  V30 BEGIN COLLECTION"))
        XCTAssertTrue(text.contains("MDLV30/HILITE ATOMS=(4 2 3 4 5) BONDS=(1 2)")
            || text.contains("MDLV30/HILITE ATOMS=(4 2 3 5 4) BONDS=(1 2)"))
        XCTAssertTrue(text.contains("M  V30 END COLLECTION"))
    }

    func testWritesStereoCollectionsAndChiralFlag() throws {
        let absGroup = CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0)
        let racGroup = CDKCxSmilesParser.encodeStereoGroup(kind: "or", number: 7)
        let relGroup = CDKCxSmilesParser.encodeStereoGroup(kind: "and", number: 5)

        let absMolecule = Molecule(
            name: "Abs",
            atoms: [
                Atom(id: 1, element: "C", position: .zero, chirality: .clockwise, cxStereoGroup: absGroup),
                Atom(id: 2, element: "O", position: CGPoint(x: 1.0, y: 0.0))
            ],
            bonds: [Bond(id: 1, a1: 1, a2: 2, order: .single)]
        )
        let absText = try CDKMDLV3000Writer.write(absMolecule)
        XCTAssertTrue(absText.contains("M  V30 COUNTS 2 1 0 0 1"))
        XCTAssertFalse(absText.contains("MDLV30/STE"))

        let mixedMolecule = Molecule(
            name: "Mixed",
            atoms: [
                Atom(id: 1, element: "C", position: .zero, chirality: .clockwise, cxStereoGroup: absGroup),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.0, y: 0.0), chirality: .anticlockwise, cxStereoGroup: racGroup),
                Atom(id: 3, element: "C", position: CGPoint(x: 2.0, y: 0.0), chirality: .clockwise, cxStereoGroup: relGroup)
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .single)
            ]
        )
        let mixedText = try CDKMDLV3000Writer.write(mixedMolecule)
        XCTAssertTrue(mixedText.contains("M  V30 COUNTS 3 2 0 0 0"))
        XCTAssertTrue(mixedText.contains("MDLV30/STEABS ATOMS=(1 1)"))
        XCTAssertTrue(mixedText.contains("MDLV30/STERAC1 ATOMS=(1 2)"))
        XCTAssertTrue(mixedText.contains("MDLV30/STEREL1 ATOMS=(1 3)"))
    }

    func testWritesDimensionFieldFor2DAnd3DScenes() throws {
        let planar = Molecule(
            name: "Planar",
            atoms: [Atom(id: 1, element: "C", position: CGPoint(x: 0.5, y: 0.5), explicitHydrogenCount: 4)],
            bonds: []
        )
        let planarText = try CDKMDLV3000Writer.write(planar)
        XCTAssertTrue(planarText.contains("2D"))

        let threeD = Molecule(
            name: "ThreeD",
            atoms: [Atom(id: 1, element: "C", position: CGPoint(x: 0.5, y: 0.5), zPosition: 0.1, explicitHydrogenCount: 4)],
            bonds: []
        )
        let threeDText = try CDKMDLV3000Writer.write(threeD)
        XCTAssertTrue(threeDText.contains("3D"))
    }

    func testWritesPositionVariationAsAttachAnyEndPoints() throws {
        let atoms = (1...9).map { index in
            Atom(id: index, element: index == 9 ? "*" : "C", position: CGPoint(x: Double(index), y: 0.0))
        }
        let bonds = [
            Bond(id: 1, a1: 1, a2: 2, order: .single),
            Bond(id: 2, a1: 2, a2: 3, order: .single),
            Bond(id: 3, a1: 3, a2: 4, order: .single),
            Bond(id: 4, a1: 4, a2: 5, order: .single),
            Bond(id: 5, a1: 5, a2: 6, order: .single),
            Bond(id: 6, a1: 6, a2: 7, order: .single),
            Bond(id: 7, a1: 7, a2: 8, order: .single),
            Bond(id: 8, a1: 8, a2: 9, order: .single)
        ]
        var molecule = Molecule(name: "PositionalVariation", atoms: atoms, bonds: bonds)
        molecule.sgroups = [
            MoleculeSgroup(kind: .extMulticenter,
                           keyword: "m",
                           atomIDs: [8, 2, 3, 4, 5, 6],
                           crossingBondIDs: [8])
        ]

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  V30 8 1 8 9 ATTACH=ANY ENDPTS=(5 2 3 4 5 6)"))
    }

    func testWritesSgroups() throws {
        let molecule = Molecule(
            name: "Sgroups",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "C", position: CGPoint(x: 1.2, y: 0.0)),
                Atom(id: 3, element: "O", position: CGPoint(x: 2.4, y: 0.0)),
                Atom(id: 4, element: "C", position: CGPoint(x: 3.6, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single),
                Bond(id: 2, a1: 2, a2: 3, order: .single),
                Bond(id: 3, a1: 3, a2: 4, order: .single)
            ]
        )

        var withSgroups = molecule
        withSgroups.sgroups = [
            MoleculeSgroup(kind: .structureRepeatUnit,
                           keyword: "SRU",
                           atomIDs: [2, 3],
                           crossingBondIDs: [1, 3],
                           subscriptText: "n",
                           roundBrackets: true,
                           connectivity: "HH",
                           brackets: [
                            MoleculeSgroupBracket(firstPoint: CGPoint(x: -1.0, y: -1.0),
                                                  secondPoint: CGPoint(x: -1.0, y: 1.0))
                           ]),
            MoleculeSgroup(kind: .generic,
                           keyword: "SUP",
                           atomIDs: [2, 3, 4],
                           crossingBondIDs: [1],
                           subscriptText: "Ph",
                           expanded: true),
            MoleculeSgroup(kind: .generic,
                           keyword: "MUL",
                           atomIDs: [2, 3],
                           crossingBondIDs: [1, 3],
                           subscriptText: "2",
                           parentAtomIDs: [2]),
            MoleculeSgroup(kind: .data,
                           keyword: "DAT",
                           atomIDs: [3],
                           crossingBondIDs: [2],
                           dataFieldName: "CDK",
                           dataValue: "AliasName")
        ]

        let text = try CDKMDLV3000Writer.write(withSgroups)

        XCTAssertTrue(text.contains("M  V30 BEGIN SGROUP"))
        XCTAssertTrue(text.contains(" SRU 0 ATOMS=(2 2 3) XBONDS=(2 1 3)"))
        XCTAssertTrue(text.contains(" LABEL=n"))
        XCTAssertTrue(text.contains(" CONNECT=HH"))
        XCTAssertTrue(text.contains(" BRKTYP=PAREN"))
        XCTAssertTrue(text.contains(" SUP 0 ATOMS=(3 2 3 4) XBONDS=(1 1) LABEL=Ph ESTATE=E"))
        XCTAssertTrue(text.contains(" MUL 0 ATOMS=(2 2 3) XBONDS=(2 1 3) MULT=2 PATOMS=(1 2)"))
        XCTAssertTrue(text.contains(" DAT 0 ATOMS=(1 3) CBONDS=(1 2) FIELDNAME=CDK FIELDDATA=AliasName"))
        XCTAssertTrue(text.contains("M  V30 END SGROUP"))
    }

    func testWritesWrappedMultipleGroupSgroups() throws {
        var atoms: [Atom] = [Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0), explicitHydrogenCount: 3)]
        for index in 2...52 {
            atoms.append(Atom(id: index, element: index == 52 ? "O" : "C", position: CGPoint(x: Double(index), y: 0.0)))
        }

        var bonds: [Bond] = []
        for index in 1...51 {
            bonds.append(Bond(id: index, a1: index, a2: index + 1, order: .single))
        }

        var molecule = Molecule(name: "MultipleGroup", atoms: atoms, bonds: bonds)
        molecule.sgroups = [
            MoleculeSgroup(kind: .generic,
                           keyword: "MUL",
                           atomIDs: Array(2...51),
                           crossingBondIDs: [1, 51],
                           subscriptText: "50",
                           parentAtomIDs: [2])
        ]

        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  V30 1 MUL 0 ATOMS=(50 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 -"))
        XCTAssertTrue(text.contains("M  V30 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 -"))
        XCTAssertTrue(text.contains("M  V30 46 47 48 49 50 51) XBONDS=(2 1 51) MULT=50 PATOMS=(1 2)"))
    }

    func testV3000WriterRoundTripsHighlightsAndTopology() throws {
        var molecule = Molecule(
            name: "RoundTrip",
            atoms: [
                Atom(id: 1, element: "C", position: CGPoint(x: 0.0, y: 0.0)),
                Atom(id: 2, element: "N", position: CGPoint(x: 1.2, y: 0.0))
            ],
            bonds: [
                Bond(id: 1, a1: 1, a2: 2, order: .single, queryType: .any, topology: .ring)
            ]
        )
        molecule.highlightedAtomIDs = [2]
        molecule.highlightedBondIDs = [1]

        let text = try CDKMDLV3000Writer.write(molecule)
        let parsed = try CDKMDLV3000Reader.read(text: text)

        XCTAssertEqual(parsed.highlightedAtomIDs, [2])
        XCTAssertEqual(parsed.highlightedBondIDs, [1])
        XCTAssertEqual(parsed.bonds.first?.queryType, .any)
        XCTAssertEqual(parsed.bonds.first?.topology, .ring)
    }

    func testRoundTripsMarkushDefinitionsFromCxSmiles() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles(
            "C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|"
        )

        let text = try CDKMDLV3000Writer.write(molecule)
        let parsed = try CDKMDLV3000Reader.read(text: text)

        XCTAssertTrue(text.contains("M  V30 BEGIN RGROUP 1"))
        XCTAssertFalse(text.contains("SUBTYPE=MARKUSHRGROUP"))
        XCTAssertEqual(parsed.atoms.filter { $0.rGroupMembership == nil && $0.symbolToDraw == "R1" }.count, 1)
        XCTAssertEqual(parsed.atoms.filter { $0.rGroupMembership == "R1" && $0.attachmentPoint == 1 }.count, 3)
        XCTAssertEqual(parsed.sgroups.count, 1)
        XCTAssertEqual(parsed.sgroups.first?.subscriptText, "1-3")

        let r1Atoms = parsed.atoms.filter { $0.rGroupMembership == "R1" }
        let grouped = Dictionary(grouping: r1Atoms, by: { $0.componentGroupID ?? -1 })

        XCTAssertEqual(grouped.count, 3)
        XCTAssertTrue(grouped.values.contains { atoms in
            atoms.contains(where: { $0.element == "O" }) && atoms.contains(where: { $0.element == "C" })
        })
        XCTAssertTrue(parsed.atoms.contains { $0.rGroupMembership == "R1" && $0.element == "Cl" })
        XCTAssertTrue(grouped.values.contains { atoms in
            atoms.contains(where: { $0.element == "C" }) && atoms.contains(where: { $0.element == "N" })
        })
    }

    func testReadsCompliantV3000RGroupBlocks() throws {
        let input = """
Molecule
  CDK     0316261200

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 R# 1.2 0 0 0 RGROUPS=(1 1)
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  V30 BEGIN RGROUP 1
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Cl 0 0 0 0 ATTCHPT=1
M  V30 END ATOM
M  V30 END CTAB
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 Br 0 0 0 0 ATTCHPT=1
M  V30 END ATOM
M  V30 END CTAB
M  V30 END RGROUP
M  END
"""

        let molecule = try CDKMDLV3000Reader.read(text: input)

        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == nil && $0.symbolToDraw == "R1" }.count, 1)
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == "R1" && $0.attachmentPoint == 1 }.count, 2)

        let grouped = Dictionary(grouping: molecule.atoms.filter { $0.rGroupMembership == "R1" },
                                 by: { $0.componentGroupID ?? -1 })
        XCTAssertEqual(grouped.count, 2)
        XCTAssertTrue(molecule.atoms.contains { $0.rGroupMembership == "R1" && $0.element == "Cl" })
        XCTAssertTrue(molecule.atoms.contains { $0.rGroupMembership == "R1" && $0.element == "Br" })
    }

    func testRoundTripsNestedSgroupHierarchyAndStereoNormalization() throws {
        let input = """

  Mrv1810 02062121432D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 11 2 0 1
M  V30 BEGIN ATOM
M  V30 1 C 0 6.16 0 0
M  V30 2 C 0 4.62 0 0 CFG=2
M  V30 3 O -1.3337 3.85 0 0
M  V30 4 C 1.3337 3.85 0 0 CFG=2
M  V30 5 O 2.6674 4.62 0 0
M  V30 6 C 1.3337 2.31 0 0
M  V30 7 C 2.6674 1.54 0 0
M  V30 8 C 2.6674 0 0 0
M  V30 9 C 1.3337 -0.77 0 0
M  V30 10 C 0 0 0 0
M  V30 11 C 0 1.54 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1
M  V30 2 1 2 3 CFG=1
M  V30 3 1 2 4
M  V30 4 1 4 5 CFG=3
M  V30 5 1 4 6
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 6 11
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 4)
M  V30 MDLV30/STEREL5 ATOMS=(1 2)
M  V30 END COLLECTION
M  V30 BEGIN SGROUP
M  V30 1 FOR 0 ATOMS=(11 1 2 3 4 5 6 7 8 9 10 11) BRKXYZ=(9 -1 -1 0 -1 8 0 0 0 0)
M  V30 2 COM 0 ATOMS=(5 1 2 3 4 5) PARENT=1 COMPNO=1 BRKXYZ=(9 2 -1 0 2 5 0 0 0 0)
M  V30 END SGROUP
M  V30 END CTAB
M  END
"""

        let molecule = try CDKMDLV3000Reader.read(text: input)
        let text = try CDKMDLV3000Writer.write(molecule)

        XCTAssertTrue(text.contains("M  V30 BEGIN COLLECTION"))
        XCTAssertTrue(text.contains("MDLV30/STEABS ATOMS=(1 4)"))
        XCTAssertTrue(text.contains("MDLV30/STEREL1 ATOMS=(1 2)"))
        XCTAssertTrue(text.contains("PARENT=1"))
        XCTAssertTrue(text.contains("COMPNO=1"))
        XCTAssertTrue(text.contains("BRKXYZ=(9"))
    }
}
