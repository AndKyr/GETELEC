import { Packet, RawData } from "engine.io-parser/build/esm/commons.js";
declare const encodePacket: ({ type, data }: Packet, supportsBinary: boolean, callback: (encodedPacket: RawData) => void) => void;
export default encodePacket;
