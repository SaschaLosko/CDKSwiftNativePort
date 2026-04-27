import Foundation
import XCTest
@testable import CDKSwiftNativePort

private struct OfficialReferenceCorpus: Decodable {
    struct Case: Decodable {
        let id: String
        let dataset: String
        let recordID: String
        let selectionReason: String
        let expectedInChI: String
        let expectedInChIKey: String
        let expectedAuxInfo: String
        let expectedExitCode: Int
        let expectedMessage: String
        let molfile: String
    }

    let schemaVersion: Int
    let suite: String
    let referenceRepository: String
    let referenceCommit: String
    let referenceSources: [String]
    let selectionLimits: [String: Int]
    let cases: [Case]
}

private struct OfficialReferenceKnownGapInventory: Codable {
    struct Gap: Codable, Hashable {
        let caseID: String
        let issue: String
    }

    let schemaVersion: Int
    let suite: String
    let referenceCommit: String
    let totalCases: Int
    let gaps: [Gap]
}

private struct OfficialReferenceGapDetail: Codable {
    let caseID: String
    let dataset: String
    let recordID: String
    let selectionReason: String
    let issue: String
    let expected: String?
    let actual: String?
    let status: String?
    let message: String?
}

private struct OfficialReferenceGapReport: Codable {
    let schemaVersion: Int
    let suite: String
    let referenceCommit: String
    let totalCases: Int
    let summary: [String: Int]
    let gaps: [OfficialReferenceGapDetail]
}

final class OfficialReferenceParityTests: XCTestCase {
    func testOfficialReferenceCorpusMetadataIsPresent() throws {
        let corpus = try loadCorpus()
        XCTAssertGreaterThanOrEqual(corpus.schemaVersion, 1)
        XCTAssertEqual(corpus.suite, "Official InChI Reference CI Corpus")
        XCTAssertEqual(corpus.referenceRepository, "https://github.com/IUPAC-InChI/InChI")
        XCTAssertEqual(corpus.referenceCommit.count, 40)
        XCTAssertEqual(corpus.selectionLimits["inchi_elemental"], 6)
        XCTAssertEqual(corpus.selectionLimits["inchi_isotopic"], 6)
        XCTAssertEqual(corpus.selectionLimits["mcule_neutral_organic"], 12)
        XCTAssertEqual(corpus.selectionLimits["mcule_charged"], 6)
        XCTAssertEqual(corpus.selectionLimits["mcule_double_stereo"], 6)
        XCTAssertEqual(corpus.selectionLimits["mcule_tetra_stereo"], 6)
        XCTAssertEqual(corpus.selectionLimits["mcule_multicomponent"], 6)
        XCTAssertEqual(corpus.cases.count, 48)
        XCTAssertEqual(Set(corpus.cases.map(\.id)).count, corpus.cases.count)
    }

    func testOfficialReferenceParityMatchesKnownGapInventory() throws {
        let corpus = try loadCorpus()
        let gaps = collectGapDetails(from: corpus)
        let report = OfficialReferenceGapReport(
            schemaVersion: 1,
            suite: corpus.suite,
            referenceCommit: corpus.referenceCommit,
            totalCases: corpus.cases.count,
            summary: Dictionary(gaps.map(\.issue).map { ($0, 1) }, uniquingKeysWith: +),
            gaps: gaps
        )
        let reportPath = try writeGapReport(report)

        if ProcessInfo.processInfo.environment["CDK_INCHI_WRITE_GAP_INVENTORY"] == "1" {
            try writeKnownGapInventory(
                OfficialReferenceKnownGapInventory(
                    schemaVersion: 1,
                    suite: corpus.suite,
                    referenceCommit: corpus.referenceCommit,
                    totalCases: corpus.cases.count,
                    gaps: gaps.map { OfficialReferenceKnownGapInventory.Gap(caseID: $0.caseID, issue: $0.issue) }
                )
            )
            return
        }

        if ProcessInfo.processInfo.environment["CDK_INCHI_REQUIRE_REFERENCE_GRADE"] == "1" {
            XCTAssertTrue(
                gaps.isEmpty,
                """
                Official InChI reference parity is not yet complete.
                Gap summary: \(formatSummary(report.summary))
                Detailed report: \(reportPath.path)
                """
            )
            return
        }

        let knownInventory = try loadKnownGapInventory()
        let currentInventory = OfficialReferenceKnownGapInventory(
            schemaVersion: 1,
            suite: corpus.suite,
            referenceCommit: corpus.referenceCommit,
            totalCases: corpus.cases.count,
            gaps: gaps.map { OfficialReferenceKnownGapInventory.Gap(caseID: $0.caseID, issue: $0.issue) }
        )

        let expectedSet = Set(knownInventory.gaps)
        let actualSet = Set(currentInventory.gaps)

        let unexpected = actualSet.subtracting(expectedSet).sorted {
            if $0.issue != $1.issue { return $0.issue < $1.issue }
            return $0.caseID < $1.caseID
        }
        let resolved = expectedSet.subtracting(actualSet).sorted {
            if $0.issue != $1.issue { return $0.issue < $1.issue }
            return $0.caseID < $1.caseID
        }

        XCTAssertTrue(
            unexpected.isEmpty && resolved.isEmpty,
            """
            Official InChI known-gap inventory drifted.
            Unexpected gaps: \(formatGapList(unexpected))
            Resolved baseline gaps: \(formatGapList(resolved))
            Current summary: \(formatSummary(report.summary))
            Detailed report: \(reportPath.path)
            Refresh inventory after review with:
            CDK_INCHI_WRITE_GAP_INVENTORY=1 swift test --package-path CDKSwiftNativePort --filter OfficialReferenceParityTests
            """
        )
    }

    private func collectGapDetails(from corpus: OfficialReferenceCorpus) -> [OfficialReferenceGapDetail] {
        var gaps: [OfficialReferenceGapDetail] = []

        for entry in corpus.cases {
            do {
                let molecule = try CDKMDLReader.read(text: entry.molfile)
                let generator = CDKInChIGeneratorFactory.shared.getInChIGenerator(molecule)
                if generator.getStatus() != .success {
                    gaps.append(
                        OfficialReferenceGapDetail(
                            caseID: entry.id,
                            dataset: entry.dataset,
                            recordID: entry.recordID,
                            selectionReason: entry.selectionReason,
                            issue: "generator_status",
                            expected: "success",
                            actual: String(describing: generator.getStatus()),
                            status: String(describing: generator.getStatus()),
                            message: generator.getMessage()
                        )
                    )
                } else {
                    let actualInChI = try generator.getInchi()
                    if actualInChI != entry.expectedInChI {
                        gaps.append(
                            OfficialReferenceGapDetail(
                                caseID: entry.id,
                                dataset: entry.dataset,
                                recordID: entry.recordID,
                                selectionReason: entry.selectionReason,
                                issue: "generator_inchi_mismatch",
                                expected: entry.expectedInChI,
                                actual: actualInChI,
                                status: String(describing: generator.getStatus()),
                                message: generator.getMessage()
                            )
                        )
                    }

                    let actualKey = try generator.getInchiKey()
                    if actualKey != entry.expectedInChIKey {
                        gaps.append(
                            OfficialReferenceGapDetail(
                                caseID: entry.id,
                                dataset: entry.dataset,
                                recordID: entry.recordID,
                                selectionReason: entry.selectionReason,
                                issue: "generator_inchikey_mismatch",
                                expected: entry.expectedInChIKey,
                                actual: actualKey,
                                status: String(describing: generator.getStatus()),
                                message: generator.getMessage()
                            )
                        )
                    }
                }
            } catch {
                gaps.append(
                    OfficialReferenceGapDetail(
                        caseID: entry.id,
                        dataset: entry.dataset,
                        recordID: entry.recordID,
                        selectionReason: entry.selectionReason,
                        issue: "generator_error",
                        expected: entry.expectedInChI,
                        actual: nil,
                        status: "error",
                        message: error.localizedDescription
                    )
                )
            }

            let parser = CDKInChIGeneratorFactory.shared.getInChIToStructure(entry.expectedInChI)
            if parser.getStatus() == .error {
                gaps.append(
                    OfficialReferenceGapDetail(
                        caseID: entry.id,
                        dataset: entry.dataset,
                        recordID: entry.recordID,
                        selectionReason: entry.selectionReason,
                        issue: "parser_error",
                        expected: entry.expectedInChI,
                        actual: nil,
                        status: String(describing: parser.getStatus()),
                        message: parser.getMessage()
                    )
                )
                continue
            }

            if parser.getStatus() == .warning {
                gaps.append(
                    OfficialReferenceGapDetail(
                        caseID: entry.id,
                        dataset: entry.dataset,
                        recordID: entry.recordID,
                        selectionReason: entry.selectionReason,
                        issue: "parser_warning",
                        expected: nil,
                        actual: nil,
                        status: String(describing: parser.getStatus()),
                        message: parser.getMessage()
                    )
                )
            }

            do {
                let reparsed = try parser.getAtomContainer()
                let roundTripGenerator = CDKInChIGeneratorFactory.shared.getInChIGenerator(reparsed)
                if roundTripGenerator.getStatus() != .success {
                    gaps.append(
                        OfficialReferenceGapDetail(
                            caseID: entry.id,
                            dataset: entry.dataset,
                            recordID: entry.recordID,
                            selectionReason: entry.selectionReason,
                            issue: "parser_roundtrip_generator_status",
                            expected: "success",
                            actual: String(describing: roundTripGenerator.getStatus()),
                            status: String(describing: roundTripGenerator.getStatus()),
                            message: roundTripGenerator.getMessage()
                        )
                    )
                } else {
                    let roundTripInChI = try roundTripGenerator.getInchi()
                    if roundTripInChI != entry.expectedInChI {
                        gaps.append(
                            OfficialReferenceGapDetail(
                                caseID: entry.id,
                                dataset: entry.dataset,
                                recordID: entry.recordID,
                                selectionReason: entry.selectionReason,
                                issue: "parser_roundtrip_inchi_mismatch",
                                expected: entry.expectedInChI,
                                actual: roundTripInChI,
                                status: String(describing: roundTripGenerator.getStatus()),
                                message: roundTripGenerator.getMessage()
                            )
                        )
                    }
                }
            } catch {
                gaps.append(
                    OfficialReferenceGapDetail(
                        caseID: entry.id,
                        dataset: entry.dataset,
                        recordID: entry.recordID,
                        selectionReason: entry.selectionReason,
                        issue: "parser_roundtrip_error",
                        expected: entry.expectedInChI,
                        actual: nil,
                        status: "error",
                        message: error.localizedDescription
                    )
                )
            }
        }

        return gaps.sorted {
            if $0.issue != $1.issue { return $0.issue < $1.issue }
            return $0.caseID < $1.caseID
        }
    }

    private func loadCorpus() throws -> OfficialReferenceCorpus {
        let url = officialReferenceDirectory().appendingPathComponent("official_reference_cases.json")
        let data = try Data(contentsOf: url)
        return try JSONDecoder().decode(OfficialReferenceCorpus.self, from: data)
    }

    private func loadKnownGapInventory() throws -> OfficialReferenceKnownGapInventory {
        let url = officialReferenceDirectory().appendingPathComponent("known_gap_inventory.json")
        let data = try Data(contentsOf: url)
        return try JSONDecoder().decode(OfficialReferenceKnownGapInventory.self, from: data)
    }

    private func writeKnownGapInventory(_ inventory: OfficialReferenceKnownGapInventory) throws {
        let url = officialReferenceDirectory().appendingPathComponent("known_gap_inventory.json")
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
        let data = try encoder.encode(inventory)
        try data.write(to: url)
    }

    private func writeGapReport(_ report: OfficialReferenceGapReport) throws -> URL {
        let url = URL(fileURLWithPath: NSTemporaryDirectory())
            .appendingPathComponent("cdk_inchi_official_reference_gap_report.json")
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
        let data = try encoder.encode(report)
        try data.write(to: url)
        return url
    }

    private func officialReferenceDirectory() -> URL {
        URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .appendingPathComponent("OfficialReference")
    }

    private func formatSummary(_ summary: [String: Int]) -> String {
        summary.keys.sorted().map { "\($0)=\(summary[$0] ?? 0)" }.joined(separator: ", ")
    }

    private func formatGapList(_ gaps: [OfficialReferenceKnownGapInventory.Gap]) -> String {
        guard !gaps.isEmpty else { return "none" }
        return gaps.prefix(12).map { "\($0.caseID):\($0.issue)" }.joined(separator: ", ")
    }
}
