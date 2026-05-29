import Foundation
import XCTest
@testable import CDKSwiftNativePort

private struct ReactionReferenceCorpus: Decodable {
    struct Case: Decodable {
        let id: String
        let format: String
        let sourceTest: String
        let sourceCase: String
        let inputFile: String
        let readerMode: String?
        let roundTrip: Bool
        let expected: ExpectedOutcome
    }

    struct ExpectedOutcome: Decodable {
        let outcome: String
        let hierarchy: ExpectedHierarchyNode?
        let errorContains: String?
    }

    struct ParticipantExpectation: Codable {
        let count: Int?
        let ids: [String]?
        let atomCounts: [Int]?
        let formulaByID: [String: String]?
    }

    struct ExpectedHierarchyNode: Codable {
        let kind: String
        let id: String?
        let isStepList: Bool?
        let reactants: ParticipantExpectation?
        let agents: ParticipantExpectation?
        let products: ParticipantExpectation?
        let properties: [String: String]?
        let children: [ExpectedHierarchyNode]?
    }

    let schemaVersion: Int
    let suite: String
    let referenceRepository: String
    let referenceTag: String
    let referenceCommit: String
    let sourceTests: [String]
    let cases: [Case]
}

private struct ReactionReferenceKnownGapInventory: Codable {
    struct Gap: Codable, Hashable {
        let caseID: String
        let phase: String
        let location: String
        let issue: String
    }

    let schemaVersion: Int
    let suite: String
    let referenceCommit: String
    let totalCases: Int
    let gaps: [Gap]
}

private struct ReactionReferenceGapDetail: Codable {
    let caseID: String
    let format: String
    let sourceTest: String
    let sourceCase: String
    let phase: String
    let location: String
    let issue: String
    let expected: String?
    let actual: String?
    let message: String?
}

private struct ReactionReferenceGapReport: Codable {
    let schemaVersion: Int
    let suite: String
    let referenceCommit: String
    let totalCases: Int
    let summary: [String: Int]
    let gaps: [ReactionReferenceGapDetail]
}

private struct ActualParticipantSnapshot: Codable {
    let ids: [String]
    let atomCounts: [Int]
    let formulasByID: [String: [String]]
}

private struct ActualHierarchySnapshot: Codable {
    let kind: String
    let id: String?
    let isStepList: Bool?
    let reactants: ActualParticipantSnapshot?
    let agents: ActualParticipantSnapshot?
    let products: ActualParticipantSnapshot?
    let properties: [String: String]?
    let children: [ActualHierarchySnapshot]
}

final class ReactionReferenceParityTests: XCTestCase {
    func testReactionReferenceCorpusMetadataIsPresent() throws {
        let corpus = try loadCorpus()
        XCTAssertGreaterThanOrEqual(corpus.schemaVersion, 1)
        XCTAssertEqual(corpus.suite, "CDK Reaction Reference Corpus")
        XCTAssertEqual(corpus.referenceRepository, "https://github.com/cdk/cdk")
        XCTAssertEqual(corpus.referenceTag, "cdk-2.12")
        XCTAssertEqual(corpus.referenceCommit.count, 40)
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("CML2Test.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("MDLRXNReaderTest.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("MDLRXNV2000ReaderTest.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("MDLRXNV3000ReaderTest.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("MDLRXNWriterTest.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("ReactionManipulatorTest.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("ReactionSetManipulatorTest.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("ReactionSchemeManipulatorTest.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("SmilesParserTest.java") }))
        XCTAssertTrue(corpus.sourceTests.contains(where: { $0.contains("CxSmilesTest.java") }))
        XCTAssertEqual(corpus.cases.count, 22)
        XCTAssertEqual(Set(corpus.cases.map(\.id)).count, corpus.cases.count)
    }

    func testReactionReferenceParityMatchesKnownGapInventory() throws {
        let corpus = try loadCorpus()
        let gaps = collectGapDetails(from: corpus)
        let report = ReactionReferenceGapReport(
            schemaVersion: 1,
            suite: corpus.suite,
            referenceCommit: corpus.referenceCommit,
            totalCases: corpus.cases.count,
            summary: Dictionary(gaps.map { ($0.issue, 1) }, uniquingKeysWith: +),
            gaps: gaps
        )
        let reportPath = try writeGapReport(report)

        if ProcessInfo.processInfo.environment["CDK_REACTION_WRITE_GAP_INVENTORY"] == "1" {
            try writeKnownGapInventory(
                ReactionReferenceKnownGapInventory(
                    schemaVersion: 1,
                    suite: corpus.suite,
                    referenceCommit: corpus.referenceCommit,
                    totalCases: corpus.cases.count,
                    gaps: gaps.map {
                        ReactionReferenceKnownGapInventory.Gap(
                            caseID: $0.caseID,
                            phase: $0.phase,
                            location: $0.location,
                            issue: $0.issue
                        )
                    }
                )
            )
            return
        }

        if ProcessInfo.processInfo.environment["CDK_REACTION_REQUIRE_UPSTREAM_PARITY"] == "1" {
            XCTAssertTrue(
                gaps.isEmpty,
                """
                CDK reaction reference parity is not yet complete.
                Gap summary: \(formatSummary(report.summary))
                Detailed report: \(reportPath.path)
                """
            )
            return
        }

        let knownInventory = try loadKnownGapInventory()
        let currentInventory = ReactionReferenceKnownGapInventory(
            schemaVersion: 1,
            suite: corpus.suite,
            referenceCommit: corpus.referenceCommit,
            totalCases: corpus.cases.count,
            gaps: gaps.map {
                ReactionReferenceKnownGapInventory.Gap(
                    caseID: $0.caseID,
                    phase: $0.phase,
                    location: $0.location,
                    issue: $0.issue
                )
            }
        )

        let expectedSet = Set(knownInventory.gaps)
        let actualSet = Set(currentInventory.gaps)
        let unexpected = actualSet.subtracting(expectedSet).sorted(by: sortGaps)
        let resolved = expectedSet.subtracting(actualSet).sorted(by: sortGaps)

        XCTAssertTrue(
            unexpected.isEmpty && resolved.isEmpty,
            """
            CDK reaction known-gap inventory drifted.
            Unexpected gaps: \(formatGapList(unexpected))
            Resolved baseline gaps: \(formatGapList(resolved))
            Current summary: \(formatSummary(report.summary))
            Detailed report: \(reportPath.path)
            Refresh inventory after review with:
            CDK_REACTION_WRITE_GAP_INVENTORY=1 swift test --package-path CDKSwiftNativePort --filter ReactionReferenceParityTests
            """
        )
    }

    private func collectGapDetails(from corpus: ReactionReferenceCorpus) -> [ReactionReferenceGapDetail] {
        var gaps: [ReactionReferenceGapDetail] = []

        for entry in corpus.cases {
            let inputURL = reactionReferenceDirectory().appendingPathComponent(entry.inputFile)
            let text: String
            do {
                text = try String(contentsOf: inputURL, encoding: .utf8)
            } catch {
                gaps.append(
                    ReactionReferenceGapDetail(
                        caseID: entry.id,
                        format: entry.format,
                        sourceTest: entry.sourceTest,
                        sourceCase: entry.sourceCase,
                        phase: "load",
                        location: "input",
                        issue: "fixture_load_error",
                        expected: entry.inputFile,
                        actual: nil,
                        message: error.localizedDescription
                    )
                )
                continue
            }

            switch entry.expected.outcome {
            case "success":
                guard let expectedHierarchy = entry.expected.hierarchy else {
                    gaps.append(
                        ReactionReferenceGapDetail(
                            caseID: entry.id,
                            format: entry.format,
                            sourceTest: entry.sourceTest,
                            sourceCase: entry.sourceCase,
                            phase: "parse",
                            location: "root",
                            issue: "missing_expected_hierarchy",
                            expected: nil,
                            actual: nil,
                            message: "Reference case declared success without an expected hierarchy."
                        )
                    )
                    continue
                }

                do {
                    let parsed = try parseReactionHierarchy(text: text, caseEntry: entry)
                    let actualSnapshot = snapshot(from: parsed)
                    gaps.append(contentsOf: compare(expected: expectedHierarchy,
                                                    actual: actualSnapshot,
                                                    caseEntry: entry,
                                                    phase: "parse",
                                                    location: "root"))

                    if entry.roundTrip {
                        do {
                            let exported = try exportReactionHierarchy(parsed, caseEntry: entry)
                            let roundTripped = try parseReactionHierarchy(text: exported, caseEntry: entry)
                            let roundTripSnapshot = snapshot(from: roundTripped)
                            gaps.append(contentsOf: compare(expected: expectedHierarchy,
                                                            actual: roundTripSnapshot,
                                                            caseEntry: entry,
                                                            phase: "roundtrip",
                                                            location: "root"))
                        } catch {
                            gaps.append(
                                ReactionReferenceGapDetail(
                                    caseID: entry.id,
                                    format: entry.format,
                                    sourceTest: entry.sourceTest,
                                    sourceCase: entry.sourceCase,
                                    phase: "roundtrip",
                                    location: "root",
                                    issue: "roundtrip_error",
                                    expected: describe(expectedHierarchy),
                                    actual: nil,
                                    message: error.localizedDescription
                                )
                            )
                        }
                    }
                } catch {
                    gaps.append(
                        ReactionReferenceGapDetail(
                            caseID: entry.id,
                            format: entry.format,
                            sourceTest: entry.sourceTest,
                            sourceCase: entry.sourceCase,
                            phase: "parse",
                            location: "root",
                            issue: "parse_error",
                            expected: describe(expectedHierarchy),
                            actual: nil,
                            message: error.localizedDescription
                        )
                    )
                }

            case "error":
                do {
                    let parsed = try parseReactionHierarchy(text: text, caseEntry: entry)
                    gaps.append(
                        ReactionReferenceGapDetail(
                            caseID: entry.id,
                            format: entry.format,
                            sourceTest: entry.sourceTest,
                            sourceCase: entry.sourceCase,
                            phase: "parse",
                            location: "root",
                            issue: "expected_error_not_thrown",
                            expected: entry.expected.errorContains,
                            actual: describe(snapshot(from: parsed)),
                            message: nil
                        )
                    )
                } catch {
                    if let expectedMessage = entry.expected.errorContains,
                       !error.localizedDescription.contains(expectedMessage) {
                        gaps.append(
                            ReactionReferenceGapDetail(
                                caseID: entry.id,
                                format: entry.format,
                                sourceTest: entry.sourceTest,
                                sourceCase: entry.sourceCase,
                                phase: "parse",
                                location: "root",
                                issue: "error_message_mismatch",
                                expected: expectedMessage,
                                actual: error.localizedDescription,
                                message: nil
                            )
                        )
                    }
                }

            default:
                gaps.append(
                    ReactionReferenceGapDetail(
                        caseID: entry.id,
                        format: entry.format,
                        sourceTest: entry.sourceTest,
                        sourceCase: entry.sourceCase,
                        phase: "parse",
                        location: "root",
                        issue: "unsupported_expected_outcome",
                        expected: entry.expected.outcome,
                        actual: nil,
                        message: nil
                    )
                )
            }
        }

        return gaps.sorted(by: sortGapDetails)
    }

    private func compare(expected: ReactionReferenceCorpus.ExpectedHierarchyNode,
                         actual: ActualHierarchySnapshot,
                         caseEntry: ReactionReferenceCorpus.Case,
                         phase: String,
                         location: String) -> [ReactionReferenceGapDetail] {
        var gaps: [ReactionReferenceGapDetail] = []

        if expected.kind != actual.kind {
            gaps.append(
                gapDetail(caseEntry: caseEntry,
                          phase: phase,
                          location: location,
                          issue: "kind_mismatch",
                          expected: expected.kind,
                          actual: actual.kind)
            )
        }

        if let expectedID = expected.id, expectedID != actual.id {
            gaps.append(
                gapDetail(caseEntry: caseEntry,
                          phase: phase,
                          location: location,
                          issue: "id_mismatch",
                          expected: expectedID,
                          actual: actual.id)
            )
        }

        if let expectedStepList = expected.isStepList, expectedStepList != actual.isStepList {
            gaps.append(
                gapDetail(caseEntry: caseEntry,
                          phase: phase,
                          location: location,
                          issue: "is_step_list_mismatch",
                          expected: String(expectedStepList),
                          actual: actual.isStepList.map { String(describing: $0) })
            )
        }

        gaps.append(contentsOf: compareParticipants(expected.reactants,
                                                    actual: actual.reactants,
                                                    role: "reactants",
                                                    caseEntry: caseEntry,
                                                    phase: phase,
                                                    location: location))
        gaps.append(contentsOf: compareParticipants(expected.agents,
                                                    actual: actual.agents,
                                                    role: "agents",
                                                    caseEntry: caseEntry,
                                                    phase: phase,
                                                    location: location))
        gaps.append(contentsOf: compareParticipants(expected.products,
                                                    actual: actual.products,
                                                    role: "products",
                                                    caseEntry: caseEntry,
                                                    phase: phase,
                                                    location: location))

        if let expectedProperties = expected.properties {
            let actualProperties = actual.properties ?? [:]
            for key in expectedProperties.keys.sorted() {
                if actualProperties[key] != expectedProperties[key] {
                    gaps.append(
                        gapDetail(caseEntry: caseEntry,
                                  phase: phase,
                                  location: location,
                                  issue: "property_mismatch",
                                  expected: "\(key)=\(expectedProperties[key] ?? "")",
                                  actual: actualProperties[key].map { "\(key)=\($0)" })
                    )
                }
            }
        }

        if let expectedChildren = expected.children {
            if expectedChildren.count != actual.children.count {
                gaps.append(
                    gapDetail(caseEntry: caseEntry,
                              phase: phase,
                              location: location,
                              issue: "children_count_mismatch",
                              expected: String(expectedChildren.count),
                              actual: String(actual.children.count))
                )
            }

            for index in 0..<min(expectedChildren.count, actual.children.count) {
                gaps.append(contentsOf: compare(expected: expectedChildren[index],
                                                actual: actual.children[index],
                                                caseEntry: caseEntry,
                                                phase: phase,
                                                location: "\(location).children[\(index)]"))
            }
        }

        return gaps
    }

    private func compareParticipants(_ expected: ReactionReferenceCorpus.ParticipantExpectation?,
                                     actual: ActualParticipantSnapshot?,
                                     role: String,
                                     caseEntry: ReactionReferenceCorpus.Case,
                                     phase: String,
                                     location: String) -> [ReactionReferenceGapDetail] {
        guard let expected else { return [] }
        guard let actual else {
            return [
                gapDetail(caseEntry: caseEntry,
                          phase: phase,
                          location: location,
                          issue: "\(role)_missing",
                          expected: describe(expected),
                          actual: nil)
            ]
        }

        var gaps: [ReactionReferenceGapDetail] = []
        if let expectedCount = expected.count, expectedCount != actual.ids.count {
            gaps.append(
                gapDetail(caseEntry: caseEntry,
                          phase: phase,
                          location: location,
                          issue: "\(role)_count_mismatch",
                          expected: String(expectedCount),
                          actual: String(actual.ids.count))
            )
        }
        if let expectedIDs = expected.ids, expectedIDs != actual.ids {
            gaps.append(
                gapDetail(caseEntry: caseEntry,
                          phase: phase,
                          location: location,
                          issue: "\(role)_ids_mismatch",
                          expected: expectedIDs.joined(separator: ","),
                          actual: actual.ids.joined(separator: ","))
            )
        }
        if let expectedAtomCounts = expected.atomCounts, expectedAtomCounts != actual.atomCounts {
            gaps.append(
                gapDetail(caseEntry: caseEntry,
                          phase: phase,
                          location: location,
                          issue: "\(role)_atom_counts_mismatch",
                          expected: expectedAtomCounts.map(String.init).joined(separator: ","),
                          actual: actual.atomCounts.map(String.init).joined(separator: ","))
            )
        }
        if let expectedFormulaByID = expected.formulaByID {
            for key in expectedFormulaByID.keys.sorted() {
                let actualFormula = actual.formulasByID[key]?.first
                if actualFormula != expectedFormulaByID[key] {
                    gaps.append(
                        gapDetail(caseEntry: caseEntry,
                                  phase: phase,
                                  location: location,
                                  issue: "\(role)_formula_mismatch",
                                  expected: "\(key)=\(expectedFormulaByID[key] ?? "")",
                                  actual: actualFormula.map { "\(key)=\($0)" })
                    )
                }
            }
        }
        return gaps
    }

    private func snapshot(from hierarchy: CDKReactionHierarchy) -> ActualHierarchySnapshot {
        switch hierarchy {
        case .reaction(let reaction):
            return snapshot(from: reaction)
        case .list(let list):
            return snapshot(fromList: list)
        case .scheme(let scheme):
            return snapshot(fromScheme: scheme)
        case .set(let set):
            return ActualHierarchySnapshot(
                kind: "set",
                id: set.id,
                isStepList: nil,
                reactants: nil,
                agents: nil,
                products: nil,
                properties: nil,
                children: set.members.map(snapshot(from:))
            )
        }
    }

    private func snapshot(from entry: CDKReactionListEntry) -> ActualHierarchySnapshot {
        switch entry {
        case .reaction(let reaction):
            return snapshot(from: reaction)
        case .list(let list):
            return snapshot(fromList: list)
        case .scheme(let scheme):
            return snapshot(fromScheme: scheme)
        }
    }

    private func snapshot(from entry: CDKReactionSchemeEntry) -> ActualHierarchySnapshot {
        switch entry {
        case .reaction(let reaction):
            return snapshot(from: reaction)
        case .list(let list):
            return snapshot(fromList: list)
        case .scheme(let scheme):
            return snapshot(fromScheme: scheme)
        }
    }

    private func snapshot(from member: CDKReactionSetMember) -> ActualHierarchySnapshot {
        switch member {
        case .reaction(let reaction):
            return snapshot(from: reaction)
        case .list(let list):
            return snapshot(fromList: list)
        case .scheme(let scheme):
            return snapshot(fromScheme: scheme)
        }
    }

    private func snapshot(fromList list: CDKReactionList) -> ActualHierarchySnapshot {
        ActualHierarchySnapshot(
            kind: "list",
            id: list.id,
            isStepList: list.isStepList,
            reactants: nil,
            agents: nil,
            products: nil,
            properties: nil,
            children: list.entries.map(snapshot(from:))
        )
    }

    private func snapshot(fromScheme scheme: CDKReactionScheme) -> ActualHierarchySnapshot {
        ActualHierarchySnapshot(
            kind: "scheme",
            id: scheme.id,
            isStepList: nil,
            reactants: nil,
            agents: nil,
            products: nil,
            properties: nil,
            children: scheme.entries.map(snapshot(from:))
        )
    }

    private func snapshot(from reaction: CDKReaction) -> ActualHierarchySnapshot {
        ActualHierarchySnapshot(
            kind: "reaction",
            id: reaction.id,
            isStepList: nil,
            reactants: participantSnapshot(reaction.reactants),
            agents: participantSnapshot(reaction.agents),
            products: participantSnapshot(reaction.products),
            properties: reaction.properties.isEmpty ? nil : reaction.properties,
            children: []
        )
    }

    private func participantSnapshot(_ molecules: [Molecule]) -> ActualParticipantSnapshot {
        var formulasByID: [String: [String]] = [:]
        let ids = molecules.map(participantLabel)
        for molecule in molecules {
            let label = participantLabel(molecule)
            let formulas = molecule.dataFieldValues(named: "Formula")
            if !formulas.isEmpty {
                formulasByID[label] = formulas
            }
        }
        return ActualParticipantSnapshot(ids: ids,
                                         atomCounts: molecules.map(\.atomCount),
                                         formulasByID: formulasByID)
    }

    private func participantLabel(_ molecule: Molecule) -> String {
        molecule.externalID ?? molecule.name
    }

    private func parseReactionHierarchy(text: String,
                                        caseEntry: ReactionReferenceCorpus.Case) throws -> CDKReactionHierarchy {
        switch caseEntry.format {
        case "rxn":
            let mode = rxnMode(for: caseEntry.readerMode)
            return hierarchy(from: try CDKRXNReader.readReactions(text: text,
                                                                  options: .init(mode: mode)))
        default:
            return try CDKFileImporter.readReactionHierarchy(text: text, fileExtension: caseEntry.format)
        }
    }

    private func exportReactionHierarchy(_ hierarchy: CDKReactionHierarchy,
                                         caseEntry: ReactionReferenceCorpus.Case) throws -> String {
        switch caseEntry.format {
        case "rxn":
            let options = CDKRXNWriter.Options(alwaysV3000: caseEntry.inputFile.contains("reaction_v3"))
            return try CDKRXNWriter.write(reactions: hierarchy.flattenedReactions,
                                          options: options)
        default:
            return try CDKFileExporter.write(reactionHierarchy: hierarchy,
                                             as: exportFormat(for: caseEntry.format))
        }
    }

    private func hierarchy(from reactions: [CDKReaction]) -> CDKReactionHierarchy {
        if reactions.count == 1, let reaction = reactions.first {
            return .reaction(reaction)
        }
        return .set(CDKReactionSet(reactions: reactions))
    }

    private func rxnMode(for rawValue: String?) -> CDKRXNReader.Mode {
        switch rawValue?.lowercased() {
        case "strict":
            return .strict
        default:
            return .relaxed
        }
    }

    private func exportFormat(for format: String) throws -> CDKFileExportFormat {
        switch format {
        case "cml":
            return .cml
        case "rxn":
            return .rxn
        case "rdf":
            return .rdf
        case "rsmi":
            return .smiles
        default:
            throw ChemError.unsupported("Unsupported reaction reference format '\(format)'.")
        }
    }

    private func loadCorpus() throws -> ReactionReferenceCorpus {
        let url = reactionReferenceDirectory().appendingPathComponent("reaction_reference_cases.json")
        let data = try Data(contentsOf: url)
        return try JSONDecoder().decode(ReactionReferenceCorpus.self, from: data)
    }

    private func loadKnownGapInventory() throws -> ReactionReferenceKnownGapInventory {
        let url = reactionReferenceDirectory().appendingPathComponent("known_gap_inventory.json")
        let data = try Data(contentsOf: url)
        return try JSONDecoder().decode(ReactionReferenceKnownGapInventory.self, from: data)
    }

    private func writeKnownGapInventory(_ inventory: ReactionReferenceKnownGapInventory) throws {
        let url = reactionReferenceDirectory().appendingPathComponent("known_gap_inventory.json")
        let data = try JSONEncoder.prettySorted.encode(inventory)
        try data.write(to: url)
    }

    private func writeGapReport(_ report: ReactionReferenceGapReport) throws -> URL {
        let data = try JSONEncoder.prettySorted.encode(report)
        let url = URL(fileURLWithPath: NSTemporaryDirectory())
            .appendingPathComponent("cdk_reaction_reference_gap_report.json")
        try data.write(to: url, options: .atomic)
        return url
    }

    private func reactionReferenceDirectory() -> URL {
        URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .appendingPathComponent("UpstreamReference")
    }

    private func gapDetail(caseEntry: ReactionReferenceCorpus.Case,
                           phase: String,
                           location: String,
                           issue: String,
                           expected: String?,
                           actual: String?,
                           message: String? = nil) -> ReactionReferenceGapDetail {
        ReactionReferenceGapDetail(
            caseID: caseEntry.id,
            format: caseEntry.format,
            sourceTest: caseEntry.sourceTest,
            sourceCase: caseEntry.sourceCase,
            phase: phase,
            location: location,
            issue: issue,
            expected: expected,
            actual: actual,
            message: message
        )
    }

    private func describe(_ node: ReactionReferenceCorpus.ExpectedHierarchyNode) -> String {
        (try? String(data: JSONEncoder.prettySorted.encode(node), encoding: .utf8)) ?? String(describing: node)
    }

    private func describe(_ participant: ReactionReferenceCorpus.ParticipantExpectation) -> String {
        (try? String(data: JSONEncoder.prettySorted.encode(participant), encoding: .utf8)) ?? String(describing: participant)
    }

    private func describe(_ node: ActualHierarchySnapshot) -> String {
        (try? String(data: JSONEncoder.prettySorted.encode(node), encoding: .utf8)) ?? String(describing: node)
    }

    private func formatSummary(_ summary: [String: Int]) -> String {
        guard !summary.isEmpty else { return "none" }
        return summary.keys.sorted().map { "\($0)=\(summary[$0] ?? 0)" }.joined(separator: ", ")
    }

    private func formatGapList(_ gaps: [ReactionReferenceKnownGapInventory.Gap]) -> String {
        guard !gaps.isEmpty else { return "none" }
        return gaps.map { "\($0.caseID):\($0.phase):\($0.location):\($0.issue)" }
            .joined(separator: ", ")
    }

    private func sortGapDetails(_ lhs: ReactionReferenceGapDetail, _ rhs: ReactionReferenceGapDetail) -> Bool {
        if lhs.caseID != rhs.caseID { return lhs.caseID < rhs.caseID }
        if lhs.phase != rhs.phase { return lhs.phase < rhs.phase }
        if lhs.location != rhs.location { return lhs.location < rhs.location }
        return lhs.issue < rhs.issue
    }

    private func sortGaps(_ lhs: ReactionReferenceKnownGapInventory.Gap, _ rhs: ReactionReferenceKnownGapInventory.Gap) -> Bool {
        if lhs.caseID != rhs.caseID { return lhs.caseID < rhs.caseID }
        if lhs.phase != rhs.phase { return lhs.phase < rhs.phase }
        if lhs.location != rhs.location { return lhs.location < rhs.location }
        return lhs.issue < rhs.issue
    }
}

private extension JSONEncoder {
    static var prettySorted: JSONEncoder {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
        return encoder
    }
}
