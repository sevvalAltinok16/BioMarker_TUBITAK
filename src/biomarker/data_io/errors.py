"""data_io hata türü. Mesajlar tek cümledir ve arayüzde olduğu gibi gösterilebilir."""


class InputError(ValueError):
    """Girdi okunamadı (dosya, kimlik veya bağlantı sorunu)."""
