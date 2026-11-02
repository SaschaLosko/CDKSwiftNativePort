import Foundation
import XCTest

final class PackageBoundaryTests: XCTestCase {
    func testSourcesContainNoAppLevelCouplingMarkers() throws {
        let root = try packageRoot()
        let sourcesRoot = root.appendingPathComponent("Sources/CDKSwiftNativePort", isDirectory: true)

        let forbiddenMarkers = [
            "import AppKit",
            "import UIKit",
            "import SwiftUI",
            "import CoreSpotlight",
            "import QuickLook",
            "import QuickLookThumbnailing",
        ]

        var violations: [String] = []
        let enumerator = FileManager.default.enumerator(
            at: sourcesRoot,
            includingPropertiesForKeys: nil)

        while let next = enumerator?.nextObject() as? URL {
            guard next.pathExtension == "swift" else { continue }
            let content = try loadUTF8Source(from: next)
            for marker in forbiddenMarkers where content.contains(marker) {
                let relativePath = next.path.replacingOccurrences(of: sourcesRoot.path + "/", with: "")
                violations.append("\(relativePath): contains '\(marker)'")
            }
        }

        XCTAssertTrue(
            violations.isEmpty,
            """
            Found package boundary violations:
            \(violations.joined(separator: "\n"))
            """)
    }

    func testPackageManifestHasNoExternalLocalPaths() throws {
        let manifest = try String(contentsOf: packageRoot().appendingPathComponent("Package.swift"), encoding: .utf8)
        XCTAssertNil(manifest.range(of: #"\.package\s*\(\s*path\s*:"#, options: .regularExpression))
        XCTAssertFalse(manifest.contains("../"), "Package targets must be self-contained.")
    }

    func testMaintainedFilesContainNoJavaSourcesOrArtifacts() throws {
        let root = try packageRoot()
        let enumerator = FileManager.default.enumerator(
            at: root, includingPropertiesForKeys: nil,
            options: [.skipsHiddenFiles])
        var violations: [String] = []
        while let file = enumerator?.nextObject() as? URL {
            if ["java", "class", "jar"].contains(file.pathExtension.lowercased()) {
                violations.append(file.path.replacingOccurrences(of: root.path + "/", with: ""))
            }
        }
        XCTAssertTrue(violations.isEmpty, "Java files are outside the standalone package contract: \(violations)")
    }

    func testPackageManifestDeclaresIOSSupport() throws {
        let root = try packageRoot()
        let manifestURL = root.appendingPathComponent("Package.swift")
        let manifest = try String(contentsOf: manifestURL, encoding: .utf8)

        XCTAssertTrue(
            manifest.contains(".iOS("),
            "Expected CDKSwiftNativePort Package.swift to declare iOS platform support.")
    }

    func testPackageManifestDeclaresDynamicLibraryProduct() throws {
        let root = try packageRoot()
        let manifestURL = root.appendingPathComponent("Package.swift")
        let manifest = try String(contentsOf: manifestURL, encoding: .utf8)

        XCTAssertTrue(
            manifest.contains("type: .dynamic"),
            "Expected CDKSwiftNativePort to remain an explicitly dynamic library product."
        )
    }

    private func packageRoot() throws -> URL {
        var url = URL(fileURLWithPath: #filePath)
        // #filePath -> .../Tests/CDKSwiftNativePortTests/PackageBoundaryTests.swift
        url.deleteLastPathComponent()
        url.deleteLastPathComponent()
        url.deleteLastPathComponent()
        guard FileManager.default.fileExists(atPath: url.appendingPathComponent("Package.swift").path) else {
            throw NSError(
                domain: "PackageBoundaryTests",
                code: 1,
                userInfo: [NSLocalizedDescriptionKey: "Could not locate package root from \(#filePath)"])
        }
        return url
    }

    private func loadUTF8Source(from url: URL) throws -> String {
        let data = try Data(contentsOf: url)
        return String(decoding: data, as: UTF8.self)
    }
}
