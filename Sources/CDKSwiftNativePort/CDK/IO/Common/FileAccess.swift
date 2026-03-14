import Foundation

public enum CDKFileAccess {
    public static func withReadableURL<T>(at url: URL,
                                          coordinateAccess: Bool = true,
                                          _ body: (URL) throws -> T) throws -> T {
        let shouldStopAccess = url.isFileURL ? url.startAccessingSecurityScopedResource() : false
        defer {
            if shouldStopAccess {
                url.stopAccessingSecurityScopedResource()
            }
        }

        guard url.isFileURL else {
            return try body(url)
        }

        if !coordinateAccess {
            return try body(url)
        }

        var coordinationError: NSError?
        var result: Result<T, Error>?
        let coordinator = NSFileCoordinator()
        coordinator.coordinate(readingItemAt: url,
                               options: [.withoutChanges],
                               error: &coordinationError) { coordinatedURL in
            result = Result {
                try body(coordinatedURL)
            }
        }

        if let coordinationError {
            throw coordinationError
        }

        guard let result else {
            throw ChemError.parseFailed("Unable to coordinate file access.")
        }
        return try result.get()
    }

    public static func readData(from url: URL,
                                coordinateAccess: Bool = true) throws -> Data {
        try withReadableURL(at: url, coordinateAccess: coordinateAccess) { readableURL in
            try Data(contentsOf: readableURL)
        }
    }

    public static func decodeText(from url: URL,
                                  preferredEncodings: [String.Encoding] = [.utf8,
                                                                           .utf16,
                                                                           .utf16LittleEndian,
                                                                           .utf16BigEndian,
                                                                           .isoLatin1],
                                  coordinateAccess: Bool = true) throws -> String {
        let data = try readData(from: url, coordinateAccess: coordinateAccess)
        for encoding in preferredEncodings {
            if let decoded = String(data: data, encoding: encoding) {
                return decoded
            }
        }
        return String(decoding: data, as: UTF8.self)
    }
}
